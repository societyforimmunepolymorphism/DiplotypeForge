"""HFX input adapter (authoritative specification).

Implements reading of the HFX (Haplotype Frequency eXchange) format as defined
by ``societyforimmunepolymorphism/hfx`` schema v0.1.1. See docs/DESIGN.md §2.

Posture (reviewer revisions 3 & 4):
  * **Schema validation is on by default** — each document is validated against
    the vendored authoritative ``hfx.schema.json`` using ``jsonschema``;
    nonconforming documents are rejected.
  * An explicit **lenient mode** (``lenient=True``) is provided for nonconforming
    historical files: schema violations and unknown keys are demoted to warnings,
    and the known typo ``outputRResolution`` is tolerated.
  * Allele identifiers are normalized to **fully qualified** form (e.g.
    ``HLA-A*02:01``) using ``nomenclatureUsed.database``.

The SIP ``hfx-tools`` package was evaluated; its validators are hand-rolled
field checks, not schema conformance, so we validate with ``jsonschema`` against
the authoritative schema instead (see DESIGN §2.2).

Locus order is taken from ``metadata.outputResolution`` and is authoritative.
"""

from __future__ import annotations

import csv
import json
import warnings
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from urllib.parse import urlparse

import jsonschema

from ..glstring import parse_haplotype
from ..model.core import Locus
from ..model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)

_KNOWN_METADATA_KEYS = frozenset(
    {
        "creationDateTime",
        "outputResolution",
        "hfeMethod",
        "cohortDescription",
        "nomenclatureUsed",
        "producer",
        "frequencyLocation",
        "frequencyFileHeader",
        "license",
        "hfxID",
    }
)

# Map HFX nomenclatureUsed.database to the GL String genetic-system prefix.
_DATABASE_PREFIX = {
    "IPD-IMGT/HLA": "HLA-",
    "IPD-KIR": "KIR",
    "IPD-MHC": "MHC-",
    "IPD-NHKIR": "NHKIR",
}


class HFXError(ValueError):
    """Raised for structurally invalid or nonconforming HFX documents."""


@dataclass(frozen=True, slots=True)
class HFXOptions:
    """Adapter options.

    ``lenient`` turns OFF the default schema-conformance enforcement, demoting
    violations to warnings for nonconforming historical files.
    """

    lenient: bool = False
    default_population_id: str = ""
    normalize_allele_ids: bool = True


def load_schema() -> dict:
    """Load the vendored authoritative HFX v0.1.1 JSON Schema."""
    text = (
        resources.files("diplotypeforge.adapters.schema")
        .joinpath("hfx.schema.json")
        .read_text(encoding="utf-8")
    )
    return json.loads(text)


def _normalized_schema() -> dict:
    """Return the authoritative schema with a known upstream defect relaxed.

    The v0.1.1 ``frequencyDataLocation`` definition is ``oneOf`` a ``uri``-format
    string and ``enum:["inline"]``. The literal ``"inline"`` validates against
    *both* branches (a plain string trivially satisfies the ``uri`` branch since
    ``format`` is an annotation, not an assertion, in JSON Schema 2020-12), so
    ``oneOf`` spuriously fails for *every* document. This is an upstream schema
    bug (it should be ``anyOf``). We relax that single ``oneOf`` to ``anyOf``
    in-memory so default conformance checking is meaningful rather than broken.
    The on-disk vendored schema is left untouched for fidelity.
    """
    schema = load_schema()
    defs = schema.get("$defs", {})
    fdl = defs.get("frequencyDataLocation")
    if isinstance(fdl, dict) and "oneOf" in fdl and "anyOf" not in fdl:
        fdl["anyOf"] = fdl.pop("oneOf")
    return schema


def _violation(lenient: bool, message: str) -> None:
    """Record a conformance problem: raise unless lenient (then warn)."""
    if lenient:
        warnings.warn(message, stacklevel=3)
    else:
        raise HFXError(message)


def validate_against_schema(doc: dict, lenient: bool) -> list[str]:
    """Validate a parsed HFX document against the authoritative schema.

    Returns a list of human-readable violation messages. In strict (default)
    mode any violation raises :class:`HFXError`; in lenient mode violations are
    collected and warned. ``format`` is treated as an annotation (not asserted),
    per JSON Schema 2020-12 defaults.
    """
    schema = _normalized_schema()
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(doc), key=lambda e: list(e.path))
    messages = []
    for err in errors:
        loc = "/".join(str(p) for p in err.path) or "<root>"
        messages.append(f"HFX schema violation at {loc}: {err.message}")
    if messages and not lenient:
        raise HFXError(
            "HFX document does not conform to authoritative schema v0.1.1:\n"
            + "\n".join(f"  {m}" for m in messages)
            + "\n(use lenient mode for nonconforming historical files)"
        )
    for m in messages:
        warnings.warn(m, stacklevel=3)
    return messages



def _resolve_loci(metadata: dict, lenient: bool) -> tuple[Locus, ...]:
    res = metadata.get("outputResolution")
    if res is None and "outputRResolution" in metadata:
        _violation(lenient, "metadata uses misspelled key 'outputRResolution'")
        res = metadata.get("outputRResolution")
    if not res:
        raise HFXError("metadata.outputResolution is required and defines locus order")
    loci = []
    for item in res:
        locus = item.get("locus")
        resolution = item.get("resolution")
        if not locus or not resolution:
            raise HFXError(f"outputResolution item missing locus/resolution: {item!r}")
        loci.append(Locus(name=locus, resolution=resolution))
    return tuple(loci)


def _check_unknown_keys(metadata: dict, lenient: bool) -> dict[str, object]:
    unknown = {k: v for k, v in metadata.items() if k not in _KNOWN_METADATA_KEYS}
    for key in unknown:
        _violation(lenient, f"unknown HFX metadata key preserved as provenance: {key!r}")
    return unknown


def _system_prefix(metadata: dict, lenient: bool) -> str:
    nomen = metadata.get("nomenclatureUsed") or {}
    database = nomen.get("database")
    prefix = _DATABASE_PREFIX.get(database)
    if prefix is None:
        _violation(
            lenient,
            f"unknown nomenclatureUsed.database {database!r}; cannot determine "
            "genetic-system prefix for fully-qualified allele ids",
        )
        return ""
    return prefix


def _qualify_token(token: str, prefix: str) -> str:
    """Prefix a bare allele token (``A*01:01g``) to fully-qualified (``HLA-A*01:01g``).

    If the token already carries a system prefix (contains ``-`` before ``*``),
    it is returned unchanged. KIR tokens (prefix ``KIR``) are concatenated
    directly since KIR locus names embed the system (e.g. ``KIR2DL1``).
    """
    if not prefix or not token:
        return token
    star = token.find("*")
    locus_part = token if star < 0 else token[:star]
    if "-" in locus_part or locus_part.startswith(prefix.rstrip("-")):
        return token
    return f"{prefix}{token}"


def _qualify_haplotype_gl(hap_gl: str, prefix: str) -> str:
    return "~".join(_qualify_token(tok, prefix) for tok in hap_gl.split("~"))


def _population_id(metadata: dict, options: HFXOptions, source_path: Path) -> tuple[str, str]:
    cohort = metadata.get("cohortDescription") or {}
    pops = cohort.get("population") or []
    if len(pops) > 1:
        _violation(
            options.lenient,
            f"HFX cohort lists {len(pops)} populations; treating file as a single "
            "population (see OPEN_QUESTIONS Q3). Using first population name as label.",
        )
    label = ""
    if pops:
        label = pops[0].get("name", "") or ""
    pop_id = options.default_population_id or label or source_path.stem
    return pop_id, label


def _qualify_loci(loci: tuple[Locus, ...], prefix: str) -> tuple[Locus, ...]:
    if not prefix:
        return loci
    out = []
    for locus in loci:
        name = locus.name
        if "-" in name or name.startswith(prefix.rstrip("-")):
            out.append(locus)
        else:
            out.append(Locus(name=f"{prefix}{name}", resolution=locus.resolution))
    return tuple(out)


def _frequency_from_inline(
    doc: dict, loci: tuple[Locus, ...], prefix: str
) -> list[FrequencyEntry]:
    entries = []
    for row in doc.get("frequencyData", []):
        hap_gl = row.get("haplotype")
        freq = row.get("frequency")
        if hap_gl is None or freq is None:
            raise HFXError(f"inline frequency row missing haplotype/frequency: {row!r}")
        hap = parse_haplotype(_qualify_haplotype_gl(hap_gl, prefix), loci)
        entries.append(FrequencyEntry(haplotype=hap, frequency=float(freq)))
    return entries


def _csv_path(frequency_location: str, base_dir: Path) -> Path:
    parsed = urlparse(frequency_location)
    if parsed.scheme in ("", "file"):
        raw = parsed.path or parsed.netloc
        if parsed.netloc and parsed.path:
            raw = parsed.netloc + parsed.path
        candidate = Path(raw)
        if not candidate.is_absolute():
            candidate = base_dir / candidate.name
        return candidate
    raise HFXError(f"unsupported frequencyLocation scheme for CSV read: {frequency_location!r}")


def _frequency_from_csv(
    metadata: dict, loci: tuple[Locus, ...], base_dir: Path, prefix: str
) -> list[FrequencyEntry]:
    location = metadata["frequencyLocation"]
    path = _csv_path(location, base_dir)
    if not path.exists():
        raise HFXError(f"frequency CSV not found: {path}")

    header_map = metadata.get("frequencyFileHeader") or {}
    hap_col = header_map.get("haplotype", "haplotype")
    freq_col = header_map.get("frequency", "frequency")
    count_col = header_map.get("count")

    entries: list[FrequencyEntry] = []
    counts_only = False
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        if hap_col not in fields:
            raise HFXError(f"CSV {path} missing haplotype column {hap_col!r}; has {fields}")
        has_freq = freq_col in fields
        has_count = bool(count_col) and count_col in fields
        if not has_freq and not has_count:
            raise HFXError(f"CSV {path} has neither frequency {freq_col!r} nor count column")
        counts_only = not has_freq and has_count

        for row in reader:
            hap = parse_haplotype(_qualify_haplotype_gl(row[hap_col], prefix), loci)
            count = float(row[count_col]) if has_count else None
            freq = float(row[freq_col]) if has_freq else (count or 0.0)
            entries.append(FrequencyEntry(haplotype=hap, frequency=freq, count=count))

    if counts_only:
        total = sum(e.count or 0.0 for e in entries) or 1.0
        warnings.warn(
            f"CSV {path} had no frequency column; deriving frequency from count "
            f"(total count {total:g})",
            stacklevel=2,
        )
        entries = [
            FrequencyEntry(haplotype=e.haplotype, frequency=(e.count or 0.0) / total, count=e.count)
            for e in entries
        ]
    return entries


def read_hfx(path: str | Path, options: HFXOptions | None = None) -> SourcePopulation:
    """Read one HFX metadata document into a :class:`SourcePopulation`.

    Validates against the authoritative schema by default; pass
    ``HFXOptions(lenient=True)`` for nonconforming historical files.
    """
    opts = options or HFXOptions()
    source_path = Path(path)
    if not source_path.exists():
        raise HFXError(f"HFX document not found: {source_path}")

    try:
        doc = json.loads(source_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise HFXError(f"invalid JSON in {source_path}: {exc}") from exc

    # Schema conformance (default on). Collected for provenance.
    schema_issues = validate_against_schema(doc, opts.lenient)

    metadata = doc.get("metadata")
    if not isinstance(metadata, dict):
        raise HFXError("HFX document missing required 'metadata' object")

    version = doc.get("version")
    if version is None:
        _violation(opts.lenient, "HFX document missing 'version' (expected '0.1.1')")
    elif version != "0.1.1":
        _violation(opts.lenient, f"HFX version {version!r} != adapter target '0.1.1'")

    loci = _resolve_loci(metadata, opts.lenient)
    unknown = _check_unknown_keys(metadata, opts.lenient)
    pop_id, label = _population_id(metadata, opts, source_path)

    prefix = _system_prefix(metadata, opts.lenient) if opts.normalize_allele_ids else ""
    loci = _qualify_loci(loci, prefix)

    location = metadata.get("frequencyLocation")
    if location is None:
        raise HFXError("metadata.frequencyLocation is required")

    if location == "inline":
        entries = _frequency_from_inline(doc, loci, prefix)
    else:
        entries = _frequency_from_csv(metadata, loci, source_path.parent, prefix)

    provenance: dict[str, object] = {
        "source_path": str(source_path),
        "hfx_version": version,
        "nomenclatureUsed": metadata.get("nomenclatureUsed"),
        "hfeMethod": metadata.get("hfeMethod"),
        "cohortDescription": metadata.get("cohortDescription"),
        "producer": metadata.get("producer"),
        "license": metadata.get("license"),
        "frequencyLocation": location,
        "allele_id_prefix": prefix,
    }
    if unknown:
        provenance["unknown_metadata"] = unknown
    if schema_issues:
        provenance["schema_issues"] = schema_issues

    dist = HaplotypeFrequencyDistribution(
        loci=loci,
        entries=tuple(entries),
        metadata=provenance,
    )
    return SourcePopulation(
        population_id=pop_id, distribution=dist, label=label, provenance=provenance
    )

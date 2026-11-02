"""Read a previously written diplotype dataset back into the model.

The standalone ``transform`` command operates on an existing ``diplotypes.tsv``
(optionally gzipped). Per-haplotype provenance is read from the explicit
``origin_1``/``origin_2`` columns — **not** reconstructed from the GL String,
which does not encode source population (OPEN_QUESTIONS Q17). The run
``manifest.json`` is read as a sidecar when available to recover the original
origin matrix and population order.
"""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

from ..glstring import parse_diplotype
from ..model.core import SimulatedIndividual


@dataclass(slots=True)
class LoadedDataset:
    individuals: list[SimulatedIndividual]
    manifest: dict | None = None


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", encoding="utf-8")
    return path.open("r", newline="", encoding="utf-8")


def read_diplotypes_tsv(path: str | Path) -> list[SimulatedIndividual]:
    """Load a diplotypes TSV into individuals, provenance from origin columns."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"diplotypes file not found: {p}")
    individuals: list[SimulatedIndividual] = []
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "diplotype_glstring", "origin_1", "origin_2"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"diplotypes file {p} missing columns: {sorted(missing)}")
        for row in reader:
            diplotype = parse_diplotype(row["diplotype_glstring"])
            individuals.append(
                SimulatedIndividual(
                    sample_id=row["sample_id"],
                    diplotype=diplotype,
                    origin_1=row["origin_1"],
                    origin_2=row["origin_2"],
                    replicate=int(row.get("replicate", 0) or 0),
                )
            )
    return individuals


def _find_manifest(diplotypes_path: Path) -> Path | None:
    candidate = diplotypes_path.parent / "manifest.json"
    return candidate if candidate.exists() else None


def load_dataset(path: str | Path) -> LoadedDataset:
    """Load diplotypes + the sibling manifest.json (if present)."""
    p = Path(path)
    individuals = read_diplotypes_tsv(p)
    manifest = None
    mpath = _find_manifest(p)
    if mpath is not None:
        manifest = json.loads(mpath.read_text(encoding="utf-8"))
    return LoadedDataset(individuals=individuals, manifest=manifest)

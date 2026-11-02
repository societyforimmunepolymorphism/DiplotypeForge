"""Phased GL String serialization and parsing.

Implements the subset of the Genotype List (GL) String grammar (Milius et al.
2013; Mack et al. 2023) needed to represent fully *phased* diplotypes, i.e.
multilocus genotypes where phase across loci is known.

Delimiter precedence (low -> high binding), per the grammar::

    '/'  allele ambiguity
    '~'  in-phase (joins alleles across loci on one chromosome => a haplotype)
    '+'  chromosome / gene-copy separator (joins the two haplotypes)
    '|'  genotype ambiguity
    '^'  locus separator

A DiplotypeForge diplotype uses only ``~`` and ``+``::

    hapA_locus1~hapA_locus2~...+hapB_locus1~hapB_locus2~...

This is the serialization contract shared with TissueTypeForge.

**Milestone-1 scope (reviewer item 4).** This module handles only the
*unambiguous phased diplotypes* DiplotypeForge produces. The ambiguity/locus
delimiters ``/`` (allele ambiguity), ``|`` (genotype ambiguity), and ``^``
(locus) are **rejected** — the full GL String grammar is intentionally not
implemented here. Allele identifiers are expected to be **fully qualified**
(genetic-system prefixed, e.g. ``HLA-A*02:01``).
"""

from __future__ import annotations

from .model.core import Allele, Diplotype, Haplotype, Locus

HAPLOTYPE_DELIM = "~"
CHROMOSOME_DELIM = "+"
# Grammar delimiters that are out of scope for Milestone 1.
UNSUPPORTED_DELIMS = ("/", "|", "^")


class GLStringError(ValueError):
    """Raised when a GL String cannot be parsed or is inconsistent with loci."""


def _reject_unsupported(gl: str) -> None:
    for delim in UNSUPPORTED_DELIMS:
        if delim in gl:
            raise GLStringError(
                f"GL String uses out-of-scope delimiter {delim!r}: {gl!r}. Milestone 1 "
                "supports only unambiguous phased diplotypes (delimiters '~' and '+')."
            )


def is_fully_qualified(token: str) -> bool:
    """Whether an allele token carries a genetic-system prefix (``SYSTEM-LOCUS*...``).

    Fully-qualified tokens contain a hyphen before the ``*`` (e.g. ``HLA-A*01:01``,
    ``KIR2DL1*0010101`` where the system is ``KIR...``). We treat the presence of a
    ``-`` in the locus portion as the qualification marker for HLA-style names.
    """
    star = token.find("*")
    locus_part = token if star < 0 else token[:star]
    return "-" in locus_part


def _locus_of_token(token: str) -> str:
    """Return the locus name portion of an allele token (before ``*``).

    ``HLA-A*01:01g`` -> ``HLA-A``;  ``HLA-DRB1*15:01`` -> ``HLA-DRB1``. A token
    without ``*`` is treated as its own locus name.
    """
    star = token.find("*")
    return token if star < 0 else token[:star]


def haplotype_to_glstring(haplotype: Haplotype) -> str:
    """Serialize a haplotype as ``allele1~allele2~...`` in locus order."""
    return HAPLOTYPE_DELIM.join(a.name for a in haplotype.alleles)


def diplotype_to_glstring(diplotype: Diplotype) -> str:
    """Serialize a phased diplotype as ``hap1+hap2`` preserving sampled order."""
    return (
        haplotype_to_glstring(diplotype.haplotype_1)
        + CHROMOSOME_DELIM
        + haplotype_to_glstring(diplotype.haplotype_2)
    )


def parse_haplotype(gl: str, loci: tuple[Locus, ...] | None = None) -> Haplotype:
    """Parse a single ``~``-joined haplotype GL String into a :class:`Haplotype`.

    If ``loci`` is provided, the allele tokens are aligned to the declared loci
    *by locus name*. Token order need not match the declared order (real HFX
    files sometimes list ``outputResolution`` in a different order than the
    haplotype tokens); alleles are reordered into the declared order so every
    haplotype in a distribution aligns positionally. A mismatch in the *set* of
    loci (missing/extra/duplicate locus) is an error. If ``loci`` is omitted,
    loci are inferred from the tokens in their given order.
    """
    gl = gl.strip()
    if not gl:
        raise GLStringError("empty haplotype GL String")
    _reject_unsupported(gl)
    if CHROMOSOME_DELIM in gl:
        raise GLStringError(
            f"haplotype must not contain {CHROMOSOME_DELIM!r}: {gl!r} (is this a diplotype?)"
        )
    tokens = gl.split(HAPLOTYPE_DELIM)
    if any(tok == "" for tok in tokens):
        raise GLStringError(f"haplotype has an empty allele token: {gl!r}")

    if loci is not None:
        if len(tokens) != len(loci):
            raise GLStringError(
                f"haplotype has {len(tokens)} loci but distribution declares "
                f"{len(loci)}: {gl!r}"
            )
        by_locus: dict[str, str] = {}
        for tok in tokens:
            tok_locus = _locus_of_token(tok)
            if tok_locus in by_locus:
                raise GLStringError(f"haplotype has duplicate locus {tok_locus!r}: {gl!r}")
            by_locus[tok_locus] = tok
        declared = [locus.name for locus in loci]
        missing = [name for name in declared if name not in by_locus]
        if missing:
            raise GLStringError(
                f"haplotype {gl!r} is missing declared loci {missing}; "
                f"has {sorted(by_locus)}, expected {declared}"
            )
        # Reorder into declared locus order so all haplotypes align positionally.
        alleles = [Allele(locus=name, name=by_locus[name]) for name in declared]
        return Haplotype(alleles=tuple(alleles))

    alleles = [Allele(locus=_locus_of_token(tok), name=tok) for tok in tokens]
    return Haplotype(alleles=tuple(alleles))


def parse_diplotype(gl: str, loci: tuple[Locus, ...] | None = None) -> Diplotype:
    """Parse a phased diplotype GL String ``hap1+hap2`` into a :class:`Diplotype`."""
    gl = gl.strip()
    if not gl:
        raise GLStringError("empty diplotype GL String")
    _reject_unsupported(gl)
    parts = gl.split(CHROMOSOME_DELIM)
    if len(parts) != 2:
        raise GLStringError(
            f"phased diplotype must have exactly two {CHROMOSOME_DELIM!r}-separated "
            f"haplotypes, found {len(parts)}: {gl!r}"
        )
    hap1 = parse_haplotype(parts[0], loci)
    hap2 = parse_haplotype(parts[1], loci)
    return Diplotype(haplotype_1=hap1, haplotype_2=hap2)

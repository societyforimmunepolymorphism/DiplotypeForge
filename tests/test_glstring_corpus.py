"""GL String round-trip corpus (interoperability contract).

A comprehensive corpus exercising the exact accepted/produced subset documented
in docs/GLSTRING_INTEROP.md. Every VALID entry must serialize->parse->serialize
to the identical string (lossless, phase- and order-preserving). Every INVALID
entry must raise GLStringError. Includes the multilocus example shared with
TissueTypeForge.
"""

from __future__ import annotations

import pytest

from diplotypeforge.glstring import (
    GLStringError,
    diplotype_to_glstring,
    haplotype_to_glstring,
    parse_diplotype,
    parse_haplotype,
)

# --------------------------------------------------------------------------- #
# Valid phased-diplotype corpus (fully qualified; the produced subset)
# --------------------------------------------------------------------------- #

VALID_DIPLOTYPES = [
    # single locus, heterozygous
    "HLA-A*01:01+HLA-A*02:01",
    # single locus, homozygous
    "HLA-A*02:01+HLA-A*02:01",
    # two loci
    "HLA-A*01:01~HLA-B*08:01+HLA-A*02:01~HLA-B*40:01",
    # the 5-locus example shared with TissueTypeForge (README output contract)
    "A*01:01g~C*07:01g~B*08:01g~DRB1*15:01~DQB1*06:02"
    "+A*02:01g~C*03:04g~B*40:01g~DRB1*04:04~DQB1*03:02g",
    # g/P/G-group and bare tokens round-trip verbatim
    "HLA-A*01:01g+HLA-A*01:01g",
    "HLA-DRB1*15:01+HLA-DRB1*04:04",
    # KIR-style tokens (no HLA hyphen; verbatim)
    "KIR2DL1*0010101~KIR2DL2*0010101+KIR2DL1*0030201~KIR2DL2*0030101",
    # fully homozygous multilocus
    "A*01:01g~B*08:01g+A*01:01g~B*08:01g",
]

VALID_HAPLOTYPES = [
    "HLA-A*02:01",
    "HLA-A*01:01~HLA-B*08:01",
    "A*01:01g~C*07:01g~B*08:01g~DRB1*15:01~DQB1*06:02",
]


@pytest.mark.parametrize("gl", VALID_DIPLOTYPES)
def test_diplotype_roundtrip_lossless(gl):
    d = parse_diplotype(gl)
    out = diplotype_to_glstring(d)
    assert out == gl  # byte-identical, phase + order preserved
    # idempotent second pass
    assert diplotype_to_glstring(parse_diplotype(out)) == gl


@pytest.mark.parametrize("gl", VALID_HAPLOTYPES)
def test_haplotype_roundtrip_lossless(gl):
    hap = parse_haplotype(gl)
    assert haplotype_to_glstring(hap) == gl


def test_phase_order_is_not_canonicalized():
    # Swapping the two haplotypes yields a DIFFERENT serialization (phase kept).
    a = "HLA-A*02:01~HLA-B*40:01"
    b = "HLA-A*01:01~HLA-B*08:01"
    d1 = parse_diplotype(f"{a}+{b}")
    d2 = parse_diplotype(f"{b}+{a}")
    assert diplotype_to_glstring(d1) == f"{a}+{b}"
    assert diplotype_to_glstring(d2) == f"{b}+{a}"
    assert diplotype_to_glstring(d1) != diplotype_to_glstring(d2)
    # but they are the same unordered genotype
    assert d1.unordered_key() == d2.unordered_key()


def test_homozygous_roundtrip_detected():
    d = parse_diplotype("A*01:01g~B*08:01g+A*01:01g~B*08:01g")
    assert d.is_homozygous
    assert diplotype_to_glstring(d) == "A*01:01g~B*08:01g+A*01:01g~B*08:01g"


# --------------------------------------------------------------------------- #
# Invalid corpus (out of subset -> GLStringError)
# --------------------------------------------------------------------------- #

INVALID = [
    "HLA-A*01:01/HLA-A*01:02",  # allele ambiguity
    "HLA-A*01:01+HLA-A*02:01|HLA-A*03:01+HLA-A*04:01",  # genotype ambiguity
    "HLA-A*01:01+HLA-A*02:01^HLA-B*08:01+HLA-B*44:01",  # locus delimiter
    "HLA-A*01:01",  # only one haplotype (no '+')
    "HLA-A*01:01+HLA-A*02:01+HLA-A*03:01",  # three haplotypes
    "",  # empty
    "HLA-A*01:01~+HLA-A*02:01",  # empty allele token
    "+HLA-A*02:01",  # empty first haplotype
]


@pytest.mark.parametrize("gl", INVALID)
def test_invalid_rejected(gl):
    with pytest.raises(GLStringError):
        parse_diplotype(gl)


def test_tissuetypeforge_shared_example_exact():
    """The exact 5-locus phased diplotype from the README output contract."""
    gl = (
        "A*01:01g~C*07:01g~B*08:01g~DRB1*15:01~DQB1*06:02"
        "+A*02:01g~C*03:04g~B*40:01g~DRB1*04:04~DQB1*03:02g"
    )
    d = parse_diplotype(gl)
    # exactly two haplotypes, 5 loci each, order preserved
    assert d.haplotype_1.n_loci == 5
    assert d.haplotype_2.n_loci == 5
    assert diplotype_to_glstring(d) == gl
    # first-locus alleles as expected (phase/order preserved)
    assert d.haplotype_1.alleles[0].name == "A*01:01g"
    assert d.haplotype_2.alleles[0].name == "A*02:01g"

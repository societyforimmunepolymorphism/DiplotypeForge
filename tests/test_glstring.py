"""GL String serialization/parsing and model invariants (acceptance test 9)."""

from __future__ import annotations

import pytest

from diplotypeforge.glstring import (
    GLStringError,
    diplotype_to_glstring,
    haplotype_to_glstring,
    is_fully_qualified,
    parse_diplotype,
    parse_haplotype,
)
from diplotypeforge.model.core import Diplotype, Locus

LOCI = (Locus("A", "g"), Locus("C", "g"), Locus("B", "g"), Locus("DRB1", "g"), Locus("DQB1", "g"))
HAP1 = "A*01:01g~C*07:01g~B*08:01g~DRB1*15:01~DQB1*06:02"
HAP2 = "A*02:01g~C*03:04g~B*40:01g~DRB1*04:04~DQB1*03:02g"


def test_haplotype_roundtrip():
    hap = parse_haplotype(HAP1, LOCI)
    assert hap.n_loci == 5
    assert haplotype_to_glstring(hap) == HAP1


def test_phased_diplotype_roundtrip_preserves_phase_and_order():
    gl = f"{HAP1}+{HAP2}"
    d = parse_diplotype(gl, LOCI)
    # order preserved
    assert haplotype_to_glstring(d.haplotype_1) == HAP1
    assert haplotype_to_glstring(d.haplotype_2) == HAP2
    # full round trip
    assert diplotype_to_glstring(d) == gl


def test_diplotype_equivalence_after_roundtrip():
    gl = f"{HAP1}+{HAP2}"
    d1 = parse_diplotype(gl, LOCI)
    d2 = parse_diplotype(diplotype_to_glstring(d1), LOCI)
    assert d1.ordered_key() == d2.ordered_key()
    assert d1.unordered_key() == d2.unordered_key()


def test_unordered_key_is_swap_invariant():
    a = parse_haplotype(HAP1, LOCI)
    b = parse_haplotype(HAP2, LOCI)
    assert Diplotype(a, b).unordered_key() == Diplotype(b, a).unordered_key()
    # but ordered keys differ
    assert Diplotype(a, b).ordered_key() != Diplotype(b, a).ordered_key()


def test_homozygous_detection():
    a = parse_haplotype(HAP1, LOCI)
    assert Diplotype(a, a).is_homozygous
    b = parse_haplotype(HAP2, LOCI)
    assert not Diplotype(a, b).is_homozygous


def test_locus_count_mismatch_raises():
    with pytest.raises(GLStringError):
        parse_haplotype("A*01:01g~B*08:01g", LOCI)  # only 2 vs 5 loci


def test_locus_name_mismatch_raises():
    bad_loci = (Locus("A", "g"), Locus("X", "g"))
    with pytest.raises(GLStringError):
        parse_haplotype("A*01:01~B*08:01", bad_loci)


def test_diplotype_must_have_two_haplotypes():
    with pytest.raises(GLStringError):
        parse_diplotype(HAP1, LOCI)  # no '+'
    with pytest.raises(GLStringError):
        parse_diplotype(f"{HAP1}+{HAP2}+{HAP1}", LOCI)


# --- Milestone-1 grammar restriction (reviewer item 4) ---------------------

@pytest.mark.parametrize("bad", [
    "HLA-A*01:01/HLA-A*01:02",            # allele ambiguity '/'
    "HLA-A*01:01+HLA-A*02:01|HLA-A*03:01+HLA-A*04:01",  # genotype ambiguity '|'
    "HLA-A*01:01+HLA-A*02:01^HLA-B*08:01+HLA-B*44:01",  # locus delimiter '^'
])
def test_rejects_out_of_scope_delimiters(bad):
    with pytest.raises(GLStringError):
        parse_diplotype(bad)


def test_fully_qualified_detection():
    assert is_fully_qualified("HLA-A*02:01")
    assert is_fully_qualified("HLA-DRB1*15:01")
    assert not is_fully_qualified("A*02:01")
    assert not is_fully_qualified("DRB1*15:01")


def test_fully_qualified_diplotype_roundtrips():
    gl = "HLA-A*01:01~HLA-B*08:01+HLA-A*02:01~HLA-B*40:01"
    d = parse_diplotype(gl)
    assert diplotype_to_glstring(d) == gl
    assert all(is_fully_qualified(a.name) for a in d.haplotype_1.alleles)

"""Frequency-distribution validation tests (revised strict defaults)."""

from __future__ import annotations

import pytest

from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
)
from diplotypeforge.validation import (
    ValidationOptions,
    aggregate_duplicates,
    renormalize,
    validate_distribution,
)
from tests.conftest import build_haplotype, build_toy_distribution


def _dist(entries):
    return HaplotypeFrequencyDistribution(loci=(Locus("A", "g"), Locus("B", "g")), entries=entries)


def test_valid_distribution_passes():
    report = validate_distribution(build_toy_distribution())
    assert report.ok
    assert report.errors == []


def test_sum_out_of_tolerance_errors_by_default():
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.5),
            FrequencyEntry(build_haplotype(("A*02:01", "B*07:02")), 0.1),
        )
    )  # sums to 0.6
    report = validate_distribution(dist)
    assert not report.ok
    assert any(i.code == "sum-out-of-tolerance" for i in report.errors)


def test_small_deviation_errors_without_authorization():
    # 0.98 total: previously auto-renormalized; now an error unless authorized.
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.5),
            FrequencyEntry(build_haplotype(("A*02:01", "B*07:02")), 0.48),
        )
    )
    assert not validate_distribution(dist).ok


def test_authorized_renormalization_warns_then_renormalizes():
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.5),
            FrequencyEntry(build_haplotype(("A*02:01", "B*07:02")), 0.48),
        )
    )
    opts = ValidationOptions(allow_renormalization=True, renorm_tolerance=0.05)
    report = validate_distribution(dist, opts)
    assert report.ok
    assert any(i.code == "sum-renormalized" for i in report.warnings)
    renormed, prov = renormalize(dist)
    assert abs(renormed.total_frequency() - 1.0) < 1e-12
    assert abs(prov.original_sum - 0.98) < 1e-12


def test_deviation_beyond_renorm_tolerance_still_errors_even_if_authorized():
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.5),
            FrequencyEntry(build_haplotype(("A*02:01", "B*07:02")), 0.1),
        )
    )  # 0.6, deviation 0.4 > renorm_tolerance
    opts = ValidationOptions(allow_renormalization=True, renorm_tolerance=0.05)
    assert not validate_distribution(dist, opts).ok


def test_negative_frequency_errors():
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 1.2),
            FrequencyEntry(build_haplotype(("A*02:01", "B*07:02")), -0.2),
        )
    )
    assert any(i.code == "negative-frequency" for i in validate_distribution(dist).errors)


def test_duplicate_haplotype_is_error_by_default():
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.5),
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.5),
        )
    )
    report = validate_distribution(dist)
    assert any(i.code == "duplicate-haplotype" for i in report.errors)


def test_aggregate_duplicates_is_explicit_and_records_provenance():
    dist = _dist(
        (
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.3),
            FrequencyEntry(build_haplotype(("A*01:01", "B*08:01")), 0.2),
            FrequencyEntry(build_haplotype(("A*02:01", "B*07:02")), 0.5),
        )
    )
    with pytest.warns(UserWarning):
        merged, prov = aggregate_duplicates(dist)
    assert merged.n_entries == 2
    assert prov.n_groups_merged == 1
    # merged frequency summed
    freqs = {e.haplotype.key(): e.frequency for e in merged.entries}
    assert abs(freqs[("A*01:01", "B*08:01")] - 0.5) < 1e-12
    # now validates clean
    assert validate_distribution(merged).ok


def test_inconsistent_loci_count_errors():
    bad_hap = Haplotype(alleles=(Allele("A", "A*01:01"),))
    dist = _dist((FrequencyEntry(bad_hap, 1.0),))
    assert any(i.code == "inconsistent-loci" for i in validate_distribution(dist).errors)

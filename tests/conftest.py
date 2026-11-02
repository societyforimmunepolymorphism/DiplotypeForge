"""Shared test fixtures and small synthetic data builders."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from diplotypeforge.model.core import Allele, Haplotype, Locus
from diplotypeforge.model.population import (
    FrequencyEntry,
    HaplotypeFrequencyDistribution,
    SourcePopulation,
)

# Two-locus toy system A, B with a handful of haplotypes summing to 1.0.
TOY_LOCI = (Locus("A", "g"), Locus("B", "g"))
TOY_FREQS = {
    ("A*01:01", "B*08:01"): 0.5,
    ("A*02:01", "B*07:02"): 0.3,
    ("A*03:01", "B*44:02"): 0.2,
}


def build_haplotype(tokens: tuple[str, ...], loci=TOY_LOCI) -> Haplotype:
    pairs = zip(tokens, loci, strict=True)
    return Haplotype(
        alleles=tuple(Allele(locus=locus.name, name=tok) for tok, locus in pairs)
    )


def build_toy_distribution() -> HaplotypeFrequencyDistribution:
    entries = tuple(
        FrequencyEntry(haplotype=build_haplotype(tokens), frequency=freq)
        for tokens, freq in TOY_FREQS.items()
    )
    return HaplotypeFrequencyDistribution(loci=TOY_LOCI, entries=entries)


def build_toy_population(pop_id: str = "TOY") -> SourcePopulation:
    return SourcePopulation(population_id=pop_id, distribution=build_toy_distribution())


@pytest.fixture
def toy_population() -> SourcePopulation:
    return build_toy_population()


@pytest.fixture
def toy_hfx_inline(tmp_path: Path) -> Path:
    """Write a tiny valid inline HFX document and return its path."""
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputResolution": [
                {"locus": "A", "resolution": "g"},
                {"locus": "B", "resolution": "g"},
            ],
            "hfeMethod": {
                "method": "EM",
                "parameters": [{"parameter": "iterations", "value": "1"}],
            },
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "TOY", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 100,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "inline",
        },
        "frequencyData": [
            {"haplotype": "A*01:01~B*08:01", "frequency": 0.5},
            {"haplotype": "A*02:01~B*07:02", "frequency": 0.3},
            {"haplotype": "A*03:01~B*44:02", "frequency": 0.2},
        ],
    }
    path = tmp_path / "toy.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


@pytest.fixture
def toy_hfx_csv(tmp_path: Path) -> Path:
    """Write a tiny CSV-mode HFX document with a custom header map + count."""
    csv_path = tmp_path / "toy.csv"
    csv_path.write_text(
        "Haplo,Count,Freq\n"
        "A*01:01~B*08:01,50,0.5\n"
        "A*02:01~B*07:02,30,0.3\n"
        "A*03:01~B*44:02,20,0.2\n",
        encoding="utf-8",
    )
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputResolution": [
                {"locus": "A", "resolution": "g"},
                {"locus": "B", "resolution": "g"},
            ],
            "hfeMethod": {
                "method": "EM",
                "parameters": [{"parameter": "iterations", "value": "1"}],
            },
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "TOY", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 100,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "file://toy.csv",
            "frequencyFileHeader": {"haplotype": "Haplo", "frequency": "Freq", "count": "Count"},
            "checkSum": "deadbeef",
        },
    }
    path = tmp_path / "toy-csv.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


@pytest.fixture
def toy_hfx_inline_unnormalized(tmp_path: Path) -> Path:
    """Inline HFX whose frequencies sum to 0.97 (needs authorized renormalization)."""
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputResolution": [
                {"locus": "A", "resolution": "g"},
                {"locus": "B", "resolution": "g"},
            ],
            "hfeMethod": {
                "method": "EM",
                "parameters": [{"parameter": "iterations", "value": "1"}],
            },
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "TOY", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 100,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "inline",
        },
        "frequencyData": [
            {"haplotype": "A*01:01~B*08:01", "frequency": 0.5},
            {"haplotype": "A*02:01~B*07:02", "frequency": 0.3},
            {"haplotype": "A*03:01~B*44:02", "frequency": 0.17},
        ],
    }
    path = tmp_path / "toy-unnorm.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path

"""HFX adapter tests (authoritative v0.1.1; schema-by-default, fully-qualified ids)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from diplotypeforge.adapters.hfx import HFXError, HFXOptions, load_schema, read_hfx


def test_load_vendored_schema():
    schema = load_schema()
    assert schema["$id"] == "https://phycus.org/hfx.schema.json"
    assert schema["properties"]["version"]["const"] == "0.1.1"


def test_read_inline_schema_valid_and_fully_qualified(toy_hfx_inline):
    pop = read_hfx(toy_hfx_inline)
    assert pop.population_id == "TOY"
    # loci normalized to fully-qualified HLA- names
    assert pop.distribution.locus_names == ("HLA-A", "HLA-B")
    assert pop.distribution.n_entries == 3
    # allele tokens fully qualified
    first = pop.distribution.entries[0].haplotype
    assert first.alleles[0].name.startswith("HLA-A*")
    assert abs(pop.distribution.total_frequency() - 1.0) < 1e-9


def test_read_csv_with_header_map_and_count(toy_hfx_csv):
    # CSV fixture carries an extra 'checkSum' metadata key -> schema-nonconforming,
    # so it requires lenient mode by default.
    pop = read_hfx(toy_hfx_csv, HFXOptions(lenient=True))
    assert pop.distribution.n_entries == 3
    counts = [e.count for e in pop.distribution.entries]
    assert set(counts) == {50.0, 30.0, 20.0}
    assert abs(pop.distribution.total_frequency() - 1.0) < 1e-9
    assert "unknown_metadata" in pop.provenance


def test_csv_with_checksum_rejected_in_strict_default(toy_hfx_csv):
    with pytest.raises(HFXError):
        read_hfx(toy_hfx_csv)  # default strict: 'checkSum' is additionalProperties


def test_locus_order_from_output_resolution(toy_hfx_inline):
    pop = read_hfx(toy_hfx_inline)
    first = pop.distribution.entries[0].haplotype
    assert first.alleles[0].locus == "HLA-A"
    assert first.alleles[1].locus == "HLA-B"


def test_typo_outputRResolution_errors_in_strict_default(tmp_path):  # noqa: N802
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputRResolution": [{"locus": "A", "resolution": "g"}],
            "hfeMethod": {"method": "EM", "parameters": []},
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "X", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 1,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "inline",
        },
        "frequencyData": [{"haplotype": "A*01:01", "frequency": 1.0}],
    }
    path = tmp_path / "typo.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(HFXError):
        read_hfx(path)  # default strict rejects the missing outputResolution


def test_typo_outputRResolution_tolerated_in_lenient(tmp_path):  # noqa: N802
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputRResolution": [
                {"locus": "A", "resolution": "g"},
                {"locus": "B", "resolution": "g"},
            ],
            "hfeMethod": {"method": "EM", "parameters": []},
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "X", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 1,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "inline",
        },
        "frequencyData": [{"haplotype": "A*01:01~B*08:01", "frequency": 1.0}],
    }
    path = tmp_path / "typo.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.warns(UserWarning):
        pop = read_hfx(path, HFXOptions(lenient=True))
    assert pop.distribution.locus_names == ("HLA-A", "HLA-B")


def test_missing_file_raises():
    with pytest.raises(HFXError):
        read_hfx("/nonexistent/path.hfx-metadata.json")


def test_official_inline_example_loads_in_strict_mode():
    # The vendored CAU example now carries 'version': '0.1.1' and conforms to the
    # authoritative schema, so it loads under default (strict) validation.
    example = Path(__file__).resolve().parents[1] / "examples" / "cau.hfx-metadata.json"
    if not example.exists():
        pytest.skip("example not present")
    pop = read_hfx(example)
    assert pop.distribution.locus_names == ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1")
    assert pop.distribution.n_entries > 50

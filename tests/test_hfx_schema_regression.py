"""Phase A regression tests for the integrated authoritative HFX schema.

These tests pin the behavior required after integrating the upstream HFX schema
correction (societyforimmunepolymorphism/hfx PR #10, merge commit
a434a217ff7c8f13a73d04b7635e34e76b3ef310):

  * the vendored schema is used VERBATIM (no in-memory oneOf->anyOf mutation);
  * valid inline HFX documents pass strict validation;
  * valid external file:// and https:// frequency *locations* pass validation;
  * invalid frequency locations are rejected;
  * missing required metadata is rejected;
  * deliberately invalid documents remain invalid;
  * historical (nonconforming) files still load via explicit lenient mode;
  * validation performs no mutation of the input document.

No network access and no dependency on the upstream repo at runtime.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import jsonschema
import pytest

from diplotypeforge.adapters.hfx import (
    HFX_SCHEMA_UPSTREAM_COMMIT,
    HFX_SCHEMA_VERSION,
    HFXError,
    HFXOptions,
    load_schema,
    read_hfx,
    schema_checksum,
    schema_provenance,
    validate_against_schema,
)


def _valid_inline_doc(location: str = "inline") -> dict:
    return {
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
            "frequencyLocation": location,
        },
        "frequencyData": [
            {"haplotype": "A*01:01~B*08:01", "frequency": 0.5},
            {"haplotype": "A*02:01~B*07:02", "frequency": 0.3},
            {"haplotype": "A*03:01~B*44:02", "frequency": 0.2},
        ],
    }


# --------------------------------------------------------------------------- #
# Schema provenance / fidelity
# --------------------------------------------------------------------------- #


def test_schema_version_unchanged_by_correction():
    schema = load_schema()
    assert schema["properties"]["version"]["const"] == "0.1.1"
    assert HFX_SCHEMA_VERSION == "0.1.1"


def test_schema_provenance_records_upstream_commit():
    prov = schema_provenance()
    assert prov["upstream_commit"] == HFX_SCHEMA_UPSTREAM_COMMIT
    assert prov["upstream_commit"] == "a434a217ff7c8f13a73d04b7635e34e76b3ef310"
    assert prov["schema_version"] == "0.1.1"
    # checksum is a 64-char hex sha256 and stable
    assert len(prov["canonical_sha256"]) == 64
    assert prov["canonical_sha256"] == schema_checksum()


def test_frequency_data_location_oneof_is_mutually_exclusive():
    """The corrected schema uses oneOf; 'inline' and a URI are mutually exclusive."""
    fdl = load_schema()["$defs"]["frequencyDataLocation"]
    assert "oneOf" in fdl and "anyOf" not in fdl
    sub = jsonschema.Draft202012Validator(fdl)
    assert sub.is_valid("inline")
    assert sub.is_valid("file://data.csv")
    assert sub.is_valid("https://example.org/freqs.csv")
    assert not sub.is_valid("not-a-uri-or-inline")


# --------------------------------------------------------------------------- #
# Valid documents pass strict validation
# --------------------------------------------------------------------------- #


def test_valid_inline_passes_strict(tmp_path: Path):
    path = tmp_path / "ok.hfx-metadata.json"
    path.write_text(json.dumps(_valid_inline_doc()), encoding="utf-8")
    pop = read_hfx(path)  # default strict
    assert pop.distribution.n_entries == 3


def test_valid_file_location_passes_strict_validation():
    # Schema-level: a file:// frequencyLocation must validate.
    doc = _valid_inline_doc(location="file://freqs.csv")
    del doc["frequencyData"]
    assert validate_against_schema(doc, lenient=False) == []


def test_valid_https_location_passes_strict_validation():
    doc = _valid_inline_doc(location="https://example.org/freqs.csv")
    del doc["frequencyData"]
    assert validate_against_schema(doc, lenient=False) == []


def test_file_location_reads_external_csv_strict(tmp_path: Path):
    (tmp_path / "freqs.csv").write_text(
        "haplotype,frequency\nA*01:01~B*08:01,0.6\nA*02:01~B*07:02,0.4\n",
        encoding="utf-8",
    )
    doc = _valid_inline_doc(location="file://freqs.csv")
    del doc["frequencyData"]
    path = tmp_path / "ext.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    pop = read_hfx(path)  # strict
    assert pop.distribution.n_entries == 2


# --------------------------------------------------------------------------- #
# Invalid documents are rejected
# --------------------------------------------------------------------------- #


def test_invalid_frequency_location_rejected(tmp_path: Path):
    doc = _valid_inline_doc(location="ftpx:/bogus")  # matches neither oneOf branch
    path = tmp_path / "bad-loc.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(HFXError):
        read_hfx(path)


def test_missing_required_metadata_rejected(tmp_path: Path):
    doc = _valid_inline_doc()
    del doc["metadata"]["nomenclatureUsed"]  # required by schema
    path = tmp_path / "missing-meta.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(HFXError):
        read_hfx(path)


def test_missing_version_rejected_strict(tmp_path: Path):
    doc = _valid_inline_doc()
    del doc["version"]
    path = tmp_path / "no-version.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(HFXError):
        read_hfx(path)


def test_additional_property_rejected_strict(tmp_path: Path):
    doc = _valid_inline_doc()
    doc["metadata"]["checkSum"] = "deadbeef"  # additionalProperties: false
    path = tmp_path / "extra.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(HFXError):
        read_hfx(path)


# --------------------------------------------------------------------------- #
# Lenient mode still works for historical files
# --------------------------------------------------------------------------- #


def test_historical_file_loads_in_lenient_mode(tmp_path: Path):
    doc = _valid_inline_doc()
    del doc["version"]  # historical pre-version file
    doc["metadata"]["checkSum"] = "deadbeef"  # legacy extra key
    path = tmp_path / "historical.hfx-metadata.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.warns(UserWarning):
        pop = read_hfx(path, HFXOptions(lenient=True))
    assert pop.distribution.n_entries == 3


# --------------------------------------------------------------------------- #
# No mutation during validation
# --------------------------------------------------------------------------- #


def test_validation_does_not_mutate_document(tmp_path: Path):
    doc = _valid_inline_doc()
    snapshot = copy.deepcopy(doc)
    validate_against_schema(doc, lenient=False)
    assert doc == snapshot


def test_load_schema_returns_fresh_unmutated_copy():
    """Repeated loads are independent; the frequencyDataLocation stays oneOf."""
    a = load_schema()
    a["$defs"]["frequencyDataLocation"]["anyOf"] = a["$defs"][
        "frequencyDataLocation"
    ].pop("oneOf")
    b = load_schema()
    assert "oneOf" in b["$defs"]["frequencyDataLocation"]
    assert "anyOf" not in b["$defs"]["frequencyDataLocation"]

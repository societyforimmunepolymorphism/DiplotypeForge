"""End-to-end Simulator/config/output tests (acceptance 9 e2e, 10)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from diplotypeforge.config import ConfigError, SimulationConfig
from diplotypeforge.glstring import parse_diplotype
from diplotypeforge.simulator import Simulator


def _write_config(tmp_path, hfx_path, **overrides):
    body = {
        "n_individuals": overrides.get("n_individuals", 2000),
        "seed": overrides.get("seed", 20120815),
        "mode": overrides.get("mode", "single_population"),
        "output_dir": str(tmp_path / "out"),
        "compress": overrides.get("compress", False),
    }
    sim_lines = "\n".join(f"{k} = {json.dumps(v)}" for k, v in body.items())
    pops = overrides.get("populations", [{"id": "TOY", "hfx": hfx_path.name}])
    pop_lines = "\n".join(
        f'[[populations]]\nid = "{p["id"]}"\nhfx = "{p["hfx"]}"' for p in pops
    )
    toml = f"[simulation]\n{sim_lines}\n\n{pop_lines}\n"
    cfg = tmp_path / "simulation.toml"
    cfg.write_text(toml, encoding="utf-8")
    return cfg


def test_end_to_end_run_and_write(tmp_path, toy_hfx_inline):
    # place config next to the hfx so relative path resolves
    cfg = _write_config(toy_hfx_inline.parent, toy_hfx_inline)
    config = SimulationConfig.from_file(cfg)
    result = Simulator(config).run()
    assert result.statistics.n_individuals == 2000
    files = result.write()
    assert all(k in files for k in ("diplotypes", "manifest", "summary"))

    # manifest records the full reproducibility set (rev 7)
    manifest = json.loads(Path(files["manifest"]).read_text())
    assert manifest["rng"]["seed"] == 20120815
    assert manifest["rng"]["generator"] == "PCG64"
    assert manifest["rng"]["replicate_seed_derivation"] == "numpy.random.SeedSequence.spawn"
    assert "software_version" in manifest
    assert "numpy_version" in manifest
    assert "python_version" in manifest
    assert "canonical_config" in manifest
    assert "reproducibility" in manifest
    # input checksums present
    assert manifest["inputs"][0]["sha256"] is not None
    assert len(manifest["inputs"][0]["sha256"]) == 64


def test_acceptance9_output_glstrings_parse_back(tmp_path, toy_hfx_inline):
    cfg = _write_config(toy_hfx_inline.parent, toy_hfx_inline, n_individuals=500)
    config = SimulationConfig.from_file(cfg)
    result = Simulator(config).run()
    files = result.write()

    # read diplotypes.tsv and parse each GL String back into a Diplotype
    lines = Path(files["diplotypes"]).read_text().strip().splitlines()
    header = lines[0].split("\t")
    gl_idx = header.index("diplotype_glstring")
    parsed = 0
    for line in lines[1:]:
        gl = line.split("\t")[gl_idx]
        d = parse_diplotype(gl)
        assert d.haplotype_1.n_loci == 2
        parsed += 1
    assert parsed == 500


def test_compressed_output(tmp_path, toy_hfx_inline):
    cfg = _write_config(toy_hfx_inline.parent, toy_hfx_inline, n_individuals=100, compress=True)
    config = SimulationConfig.from_file(cfg)
    result = Simulator(config).run()
    files = result.write()
    assert files["diplotypes"].endswith(".gz")
    with gzip.open(files["diplotypes"], "rt") as fh:
        assert fh.readline().startswith("sample_id")


def test_reproducible_end_to_end(tmp_path, toy_hfx_inline):
    cfg = _write_config(toy_hfx_inline.parent, toy_hfx_inline, n_individuals=1000)
    config = SimulationConfig.from_file(cfg)
    r1 = Simulator(config).run()
    r2 = Simulator(config).run()
    assert [i.sample_id for i in r1.individuals] == [i.sample_id for i in r2.individuals]
    assert [i.diplotype.ordered_key() for i in r1.individuals] == [
        i.diplotype.ordered_key() for i in r2.individuals
    ]


# --- acceptance 10: invalid configurations fail with informative errors ----

def test_acceptance10_unknown_mode_errors():
    with pytest.raises(ConfigError, match="unknown mode"):
        SimulationConfig.from_dict(
            {"simulation": {"mode": "bogus"}, "populations": [{"id": "X", "hfx": "x.json"}]}
        )


def test_acceptance10_no_populations_errors():
    with pytest.raises(ConfigError):
        SimulationConfig.from_dict({"simulation": {"n_individuals": 10}})


def test_acceptance10_single_pop_requires_one_population():
    with pytest.raises(ConfigError, match="exactly one population"):
        SimulationConfig.from_dict(
            {
                "simulation": {"mode": "single_population"},
                "populations": [{"id": "A", "hfx": "a.json"}, {"id": "B", "hfx": "b.json"}],
            }
        )


def test_acceptance10_wahlund_requires_two_populations():
    with pytest.raises(ConfigError, match="wahlund mode requires"):
        SimulationConfig.from_dict(
            {"simulation": {"mode": "wahlund"}, "populations": [{"id": "X", "hfx": "x.json"}]}
        )


def test_acceptance10_admixture_requires_a():
    with pytest.raises(ConfigError, match="admixture_a"):
        SimulationConfig.from_dict(
            {
                "simulation": {"mode": "admixture"},
                "populations": [{"id": "A", "hfx": "a.json"}, {"id": "B", "hfx": "b.json"}],
            }
        )


def test_acceptance10_matrix_dims_must_match_populations():
    with pytest.raises(ConfigError, match="dims must match"):
        SimulationConfig.from_dict(
            {
                "simulation": {"mode": "matrix", "matrix": [[1.0]]},
                "populations": [{"id": "A", "hfx": "a.json"}, {"id": "B", "hfx": "b.json"}],
            }
        )


def test_acceptance10_negative_n_errors():
    with pytest.raises(ConfigError):
        SimulationConfig.from_dict(
            {
                "simulation": {"n_individuals": -5, "mode": "single_population"},
                "populations": [{"id": "X", "hfx": "x.json"}],
            }
        )


def test_authorized_renormalization_recorded_in_manifest(tmp_path, toy_hfx_inline_unnormalized):
    cfg = _write_config(
        toy_hfx_inline_unnormalized.parent,
        toy_hfx_inline_unnormalized,
        n_individuals=200,
    )
    # inject authorization into the TOML
    text = cfg.read_text() + '\n'
    text = text.replace(
        "[simulation]",
        "[simulation]\nallow_renormalization = true\nrenorm_tolerance = 0.05",
    )
    cfg.write_text(text)
    config = SimulationConfig.from_file(cfg)
    result = Simulator(config).run()
    files = result.write()
    manifest = json.loads(Path(files["manifest"]).read_text())
    assert "renormalization" in manifest
    rec = next(iter(manifest["renormalization"].values()))
    assert rec["original_sum"] < 1.0
    assert abs(rec["scale_factor"] * rec["original_sum"] - 1.0) < 1e-9

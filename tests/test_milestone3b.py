"""Milestone 3B tests: reproducible replicates and aggregate statistics.

Covers brief §10 items 11, 12, 13 plus the aggregate-statistics contract.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from diplotypeforge.config import SimulationConfig, TransformationConfig
from diplotypeforge.replicates import aggregate_replicates
from diplotypeforge.simulator import Simulator


def _write_hfx(tmp_path: Path) -> Path:
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputResolution": [{"locus": "A", "resolution": "g"}],
            "hfeMethod": {"method": "direct", "parameters": []},
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "POP", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 100,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "inline",
        },
        "frequencyData": [
            {"haplotype": "A*01:01", "frequency": 0.4},
            {"haplotype": "A*02:01", "frequency": 0.3},
            {"haplotype": "A*03:01", "frequency": 0.2},
            {"haplotype": "A*11:01", "frequency": 0.1},
        ],
    }
    p = tmp_path / "pop.hfx-metadata.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def _config(tmp_path, **over) -> SimulationConfig:
    hfx = _write_hfx(tmp_path)
    data = {
        "simulation": {
            "n_individuals": over.pop("n", 1000),
            "seed": over.pop("seed", 20120815),
            "mode": "single_population",
            "replicates": over.pop("replicates", 5),
            "output_dir": str(tmp_path / "out"),
        },
        "populations": [{"id": "POP", "hfx": hfx.name}],
    }
    if "transformation" in over:
        data["transformation"] = over.pop("transformation")
    return SimulationConfig.from_dict(data, base_dir=str(tmp_path))


def test_acceptance11_replicates_independent_and_reproducible(tmp_path):
    cfg = _config(tmp_path, replicates=5)
    r1 = Simulator(cfg).run()
    r2 = Simulator(cfg).run()
    # Reproducible: same keys across the whole batch.
    k1 = [i.diplotype.ordered_key() for i in r1.individuals]
    k2 = [i.diplotype.ordered_key() for i in r2.individuals]
    assert k1 == k2
    # Independent: replicate 0 and replicate 1 are different draws.
    n = cfg.n_individuals
    rep0 = [i.diplotype.ordered_key() for i in r1.individuals[:n]]
    rep1 = [i.diplotype.ordered_key() for i in r1.individuals[n : 2 * n]]
    assert rep0 != rep1


def test_acceptance12_single_replicate_matches_batch(tmp_path):
    # Replicate index 4 ("replicate 5") must be identical whether generated
    # alone or inside a batch of 100.
    cfg_small = _config(tmp_path, replicates=5)
    cfg_big = replace(cfg_small, replicates=100)

    sim_small = Simulator(cfg_small)
    sim_small.validate()
    from diplotypeforge.sampling.origin import OriginModel
    from diplotypeforge.sampling.sampler import Sampler

    origin = OriginModel.single_population()
    sampler = Sampler(sim_small.source_populations, origin)
    # replicate 4 generated directly for each config's Simulator
    _, final_small, _ = sim_small.generate_replicate(4, origin, sampler, "DF-x")

    sim_big = Simulator(cfg_big)
    sim_big.validate()
    sampler_big = Sampler(sim_big.source_populations, origin)
    _, final_big, _ = sim_big.generate_replicate(4, origin, sampler_big, "DF-x")

    assert [i.diplotype.ordered_key() for i in final_small] == [
        i.diplotype.ordered_key() for i in final_big
    ]


def test_acceptance13_sampled_unchanged_transform_on_off(tmp_path):
    base = _config(tmp_path, replicates=3)
    on = replace(
        base,
        transformation=TransformationConfig(
            enabled=True, target_fraction=0.5, strategy="reassignment"
        ),
    )
    off = replace(base, transformation=TransformationConfig(enabled=False))
    r_on = Simulator(on).run()
    r_off = Simulator(off).run()
    # The sampled (pre-transform) diplotypes must be identical regardless of
    # whether the transform is enabled (separate RNG streams).
    on_sampled = [i.diplotype.ordered_key() for i in r_on.original_individuals]
    off_sampled = [i.diplotype.ordered_key() for i in r_off.individuals]
    assert on_sampled == off_sampled
    # And the transform did change the final sample.
    on_final = [i.diplotype.ordered_key() for i in r_on.individuals]
    assert on_final != on_sampled


def test_aggregate_statistics_shape_and_counts(tmp_path):
    cfg = _config(
        tmp_path,
        replicates=6,
        transformation={
            "enabled": True,
            "type": "homozygosity",
            "strategy": "reassignment",
            "target_fraction": 0.4,
        },
    )
    result = Simulator(cfg).run()
    bs = result.replicate_summary
    assert bs is not None
    assert bs.n_replicates == 6
    # Each enabled replicate that reached its target counts as success.
    assert bs.n_success + bs.n_fail == 6
    # Aggregates contain mean/sd/median/quantiles for the key fields.
    agg = bs.aggregates["observed_homozygosity"]
    for k in ("mean", "sd", "median", "q05", "q25", "q75", "q95"):
        assert k in agg
    # Machine-readable round-trip.
    payload = bs.to_dict()
    assert json.loads(json.dumps(payload))["n_replicates"] == 6
    assert len(payload["per_replicate"]) == 6


def test_replicates_json_written(tmp_path):
    cfg = _config(tmp_path, replicates=4)
    result = Simulator(cfg).run()
    files = result.write()
    assert "replicates" in files
    data = json.loads(Path(files["replicates"]).read_text())
    assert data["n_replicates"] == 4
    assert "aggregates" in data


def test_aggregate_empty_is_safe():
    bs = aggregate_replicates([])
    assert bs.n_replicates == 0
    assert bs.n_success == 0

"""Phase B: provenance and reproducibility verification.

Goes beyond matching summary statistics: reproducibility is asserted via
**byte-level checksums of the written output files** (diplotypes.tsv,
*_frequencies.tsv, summary.json), and provenance is asserted by tracking each
haplotype copy's source-population origin through the transform.

Covered:
  * source-population identity travels with each haplotype copy on reassignment;
  * original vs transformed origin-pair distributions are recorded separately;
  * the transformed sample is never described as sampled from M;
  * standalone transform preserves/report provenance;
  * input HFX checksums + schema provenance are recorded in the manifest;
  * identical config/inputs/seed -> byte-identical output checksums;
  * replicate k alone == replicate k in a batch (checksum of its slice);
  * enabling the transform does not change the sampled diplotypes;
  * serialization is deterministic;
  * input ordering (population order) does not change a single-population run.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np

from diplotypeforge.config import SimulationConfig, TransformationConfig
from diplotypeforge.model.core import Allele, Diplotype, Haplotype, SimulatedIndividual
from diplotypeforge.simulator import Simulator
from diplotypeforge.transform import TransformationTarget, transform_homozygosity
from diplotypeforge.transform.diagnostics import transform_diagnostics
from diplotypeforge.transform.feasibility import haplotype_counts


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()


def _write_hfx(path: Path, name: str, freqs: dict[str, float]) -> Path:
    doc = {
        "version": "0.1.1",
        "metadata": {
            "outputResolution": [{"locus": "A", "resolution": "g"}],
            "hfeMethod": {"method": "direct", "parameters": []},
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": name, "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 100,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.20.0"},
            "frequencyLocation": "inline",
        },
        "frequencyData": [{"haplotype": h, "frequency": f} for h, f in freqs.items()],
    }
    p = path / f"{name}.hfx-metadata.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def _single_config(tmp_path: Path, **over) -> SimulationConfig:
    _write_hfx(tmp_path, "POP", {"A*01:01": 0.4, "A*02:01": 0.3, "A*03:01": 0.2, "A*11:01": 0.1})
    data = {
        "simulation": {
            "n_individuals": over.pop("n", 500),
            "seed": over.pop("seed", 20120815),
            "mode": "single_population",
            "replicates": over.pop("replicates", 3),
            "output_dir": str(tmp_path / over.pop("outdir", "out")),
        },
        "populations": [{"id": "POP", "hfx": "POP.hfx-metadata.json"}],
    }
    if "transformation" in over:
        data["transformation"] = over.pop("transformation")
    return SimulationConfig.from_dict(data, base_dir=str(tmp_path))


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #


def test_origin_travels_with_haplotype_copies_through_transform():
    """Each distinct source-origin label count is conserved by reassignment."""
    inds = [
        SimulatedIndividual(
            sample_id=f"s{i}",
            diplotype=Diplotype(
                haplotype_1=Haplotype((Allele("A", f"A*{i % 4:02d}:01"),)),
                haplotype_2=Haplotype((Allele("A", f"A*{(i + 1) % 4:02d}:01"),)),
            ),
            origin_1="P0" if i % 2 == 0 else "P1",
            origin_2="P1" if i % 3 == 0 else "P0",
        )
        for i in range(40)
    ]

    def origin_multiset(pop):
        c = Counter()
        for x in pop:
            c[x.origin_1] += 1
            c[x.origin_2] += 1
        return c

    def key_origin_pairs(pop):
        # (haplotype key, origin) multiset must be conserved: origin rides the copy.
        c = Counter()
        for x in pop:
            c[(x.diplotype.haplotype_1.key(), x.origin_1)] += 1
            c[(x.diplotype.haplotype_2.key(), x.origin_2)] += 1
        return c

    before_origins = origin_multiset(inds)
    before_pairs = key_origin_pairs(inds)
    transformed, _ = transform_homozygosity(
        inds, TransformationTarget(fraction=0.5), rng=np.random.default_rng(1)
    )
    assert origin_multiset(transformed) == before_origins
    # The (key, origin) multiset is conserved: the origin label stays attached
    # to the exact haplotype copy it came with.
    assert key_origin_pairs(transformed) == before_pairs


def test_original_and_transformed_pair_distributions_recorded_separately():
    inds = [
        SimulatedIndividual(
            sample_id=f"s{i}",
            diplotype=Diplotype(
                haplotype_1=Haplotype((Allele("A", f"A*{i % 3:02d}:01"),)),
                haplotype_2=Haplotype((Allele("A", f"A*{(i + 1) % 3:02d}:01"),)),
            ),
            origin_1="P0",
            origin_2="P1",
        )
        for i in range(30)
    ]
    transformed, _ = transform_homozygosity(
        inds, TransformationTarget(fraction=0.5), rng=np.random.default_rng(2)
    )
    # single-locus trivial pops for the diagnostics call
    from diplotypeforge.model.core import Locus
    from diplotypeforge.model.population import (
        FrequencyEntry,
        HaplotypeFrequencyDistribution,
        SourcePopulation,
    )

    dist = HaplotypeFrequencyDistribution(
        loci=(Locus("A", "g"),),
        entries=(
            FrequencyEntry(Haplotype((Allele("A", "A*00:01"),)), 0.34),
            FrequencyEntry(Haplotype((Allele("A", "A*01:01"),)), 0.33),
            FrequencyEntry(Haplotype((Allele("A", "A*02:01"),)), 0.33),
        ),
    )
    pops = [SourcePopulation("P0", dist), SourcePopulation("P1", dist)]
    diag = transform_diagnostics(inds, transformed, pops, np.array([[0.5, 0.0], [0.0, 0.5]]))
    # before/after origin-pair distributions are separate fields
    assert diag.origin_pair_before is not diag.origin_pair_after
    # and the diagnostics explicitly state the transformed sample is NOT ~ M
    assert "NOT distributed according to the original" in diag.note


def test_manifest_records_input_checksums_and_schema_provenance(tmp_path):
    cfg = _single_config(tmp_path, replicates=1)
    result = Simulator(cfg).run()
    m = result.manifest
    assert m["inputs"][0]["sha256"] and len(m["inputs"][0]["sha256"]) == 64
    assert m["hfx_schema"]["upstream_commit"] == "a434a217ff7c8f13a73d04b7635e34e76b3ef310"
    assert len(m["hfx_schema"]["canonical_sha256"]) == 64


# --------------------------------------------------------------------------- #
# Reproducibility by OUTPUT CHECKSUM
# --------------------------------------------------------------------------- #

_OUTPUT_FILES = (
    "diplotypes.tsv",
    "haplotype_frequencies.tsv",
    "diplotype_frequencies.tsv",
    "summary.json",
)


def _run_and_checksum(cfg: SimulationConfig, out: Path) -> dict[str, str]:
    result = Simulator(cfg).run()
    result.write(out)
    return {name: _sha256_file(out / name) for name in _OUTPUT_FILES}


def test_identical_config_gives_identical_output_checksums(tmp_path):
    cfg = _single_config(tmp_path, replicates=3)
    a = _run_and_checksum(cfg, tmp_path / "a")
    b = _run_and_checksum(cfg, tmp_path / "b")
    assert a == b


def test_transform_run_is_reproducible_by_checksum(tmp_path):
    cfg = _single_config(
        tmp_path,
        replicates=3,
        transformation={
            "enabled": True,
            "type": "homozygosity",
            "strategy": "reassignment",
            "target_fraction": 0.5,
        },
    )
    a = _run_and_checksum(cfg, tmp_path / "ta")
    b = _run_and_checksum(cfg, tmp_path / "tb")
    assert a == b


def test_enabling_transform_does_not_change_sampled_diplotypes_checksum(tmp_path):
    base = _single_config(tmp_path, replicates=3)
    off = replace(base, transformation=TransformationConfig(enabled=False))
    on = replace(
        base,
        transformation=TransformationConfig(
            enabled=True, target_fraction=0.5, strategy="reassignment"
        ),
    )
    # off run's diplotypes.tsv == on run's diplotypes_original.tsv
    r_off = Simulator(off).run()
    r_off.write(tmp_path / "off")
    r_on = Simulator(on).run()
    files_on = r_on.write(tmp_path / "on")
    off_sum = _sha256_file(tmp_path / "off" / "diplotypes.tsv")
    on_original_sum = _sha256_file(Path(files_on["diplotypes_original"]))
    assert off_sum == on_original_sum
    # and the transformed output differs from the original
    assert _sha256_file(tmp_path / "on" / "diplotypes.tsv") != on_original_sum


def test_replicate_slice_checksum_alone_equals_in_batch(tmp_path):
    """Replicate k's diplotype draw is identical whether k is alone or in a batch.

    Note: the per-individual ``sample_id`` embeds the run hash, which binds the
    full canonical config (including the replicate count), so sample_ids legitimately
    differ between a 1-replicate and a 10-replicate run. Reproducibility of the
    *scientific content* is therefore asserted on the diplotype + origin columns,
    serialized without the run-hash-dependent sample_id.
    """
    small = _single_config(tmp_path, replicates=1, seed=42, outdir="small")
    big = replace(small, replicates=10, output_dir=str(tmp_path / "big"))

    sim_small = Simulator(small)
    r_small = sim_small.run()
    n = small.n_individuals
    rep0_small = r_small.individuals[:n]

    sim_big = Simulator(big)
    r_big = sim_big.run()
    rep0_big = r_big.individuals[:n]

    def content_checksum(pop):
        h = hashlib.sha256()
        for ind in pop:
            row = (
                f"{ind.replicate}\t"
                f"{ind.diplotype.ordered_key()}\t{ind.origin_1}\t{ind.origin_2}\n"
            )
            h.update(row.encode("utf-8"))
        return h.hexdigest()

    assert content_checksum(rep0_small) == content_checksum(rep0_big)


def test_serialization_is_deterministic(tmp_path):
    cfg = _single_config(tmp_path, replicates=2)
    result = Simulator(cfg).run()
    from diplotypeforge.output import write_diplotypes

    p1 = tmp_path / "d1.tsv"
    p2 = tmp_path / "d2.tsv"
    write_diplotypes(result.individuals, p1)
    write_diplotypes(result.individuals, p2)
    assert _sha256_file(p1) == _sha256_file(p2)


def test_input_population_label_does_not_change_single_population_output(tmp_path):
    """Re-running after rewriting the HFX with reordered frequencyData rows must
    not change outputs (haplotypes are sorted into a frozen array)."""
    # original order
    _write_hfx(
        tmp_path, "ORD", {"A*01:01": 0.4, "A*02:01": 0.3, "A*03:01": 0.2, "A*11:01": 0.1}
    )
    data = {
        "simulation": {
            "n_individuals": 500,
            "seed": 7,
            "mode": "single_population",
            "replicates": 2,
            "output_dir": str(tmp_path / "o1"),
        },
        "populations": [{"id": "ORD", "hfx": "ORD.hfx-metadata.json"}],
    }
    cfg1 = SimulationConfig.from_dict(data, base_dir=str(tmp_path))
    a = _run_and_checksum(cfg1, tmp_path / "o1")

    # reordered frequencyData (same multiset of entries)
    _write_hfx(
        tmp_path, "ORD", {"A*11:01": 0.1, "A*03:01": 0.2, "A*01:01": 0.4, "A*02:01": 0.3}
    )
    data2 = dict(data)
    data2["simulation"] = dict(data["simulation"], output_dir=str(tmp_path / "o2"))
    cfg2 = SimulationConfig.from_dict(data2, base_dir=str(tmp_path))
    b = _run_and_checksum(cfg2, tmp_path / "o2")
    assert a == b


def test_counts_preserved_end_to_end_with_transform(tmp_path):
    cfg = _single_config(
        tmp_path,
        replicates=1,
        n=2000,
        transformation={
            "enabled": True,
            "type": "homozygosity",
            "strategy": "reassignment",
            "target_fraction": 0.3,
        },
    )
    result = Simulator(cfg).run()
    before = Counter(dict(haplotype_counts(result.original_individuals)))
    after = Counter(dict(haplotype_counts(result.individuals)))
    assert before == after

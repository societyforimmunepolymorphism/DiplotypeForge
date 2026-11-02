"""High-level Simulator and SimulationResult (the primary Python API)."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version
from pathlib import Path

import numpy as np

from .adapters.hfx import HFXOptions, read_hfx, schema_provenance
from .config import SimulationConfig
from .expectation import marginal_ancestry_proportions
from .model.core import SimulatedIndividual
from .model.population import SourcePopulation
from .output import (
    write_diplotype_frequencies,
    write_diplotypes,
    write_haplotype_frequencies,
    write_manifest,
    write_summary,
)
from .replicates import (
    ReplicateBatchSummary,
    ReplicateStatistics,
    aggregate_replicates,
)
from .sampling.origin import OriginModel
from .sampling.sampler import Sampler, make_rng
from .stats import SummaryStatistics, compute_summary
from .transform import (
    TransformationResult,
    TransformationTarget,
    transform_diagnostics,
    transform_homozygosity,
)
from .validation import (
    AggregationProvenance,
    Renormalization,
    ValidationOptions,
    aggregate_duplicates,
    renormalize,
    validate_distribution,
)


def _software_version() -> str:
    try:
        return _pkg_version("diplotypeforge")
    except PackageNotFoundError:
        return "0.1.0"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _canonical_config(config: SimulationConfig) -> dict:
    """A normalized, order-stable view of the configuration for hashing + manifest."""
    return {
        "n_individuals": config.n_individuals,
        "seed": config.seed,
        "mode": config.mode,
        "replicates": config.replicates,
        "compress": config.compress,
        "strict": config.strict,
        "lenient_hfx": config.lenient_hfx,
        "allow_renormalization": config.allow_renormalization,
        "renorm_tolerance": config.renorm_tolerance,
        "aggregate_duplicates": config.aggregate_duplicates,
        "populations": [{"id": p.id, "hfx": p.hfx} for p in config.populations],
    }


def _run_hash(canonical_config: dict, software_version: str) -> str:
    payload = {"config": canonical_config, "software_version": software_version}
    blob = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:8]


def _build_origin_model(config: SimulationConfig, n_pops: int) -> OriginModel:
    mode = config.mode
    if mode == "single_population":
        return OriginModel.single_population()
    if mode == "wahlund":
        weights = list(config.weights) if config.weights is not None else [1.0] * n_pops
        return OriginModel.wahlund(weights)
    if mode == "admixture":
        return OriginModel.admixture(float(config.admixture_a))
    if mode == "matrix":
        return OriginModel.from_matrix(
            [list(row) for row in config.matrix], symmetric=config.symmetric
        )
    raise NotImplementedError(f"unsupported mode {mode!r}")


def _haplotype_freq_rmse(stats: SummaryStatistics) -> float:
    """RMS of (observed - input) haplotype frequency over the comparison set."""
    comps = stats.haplotype_comparisons
    if not comps:
        return 0.0
    sq = sum((c.observed_frequency - c.input_frequency) ** 2 for c in comps)
    return (sq / len(comps)) ** 0.5


def _replicate_record(
    rep: int,
    stats: SummaryStatistics,
    transform_enabled: bool,
    tresult,
) -> ReplicateStatistics:
    return ReplicateStatistics(
        replicate=rep,
        n_individuals=stats.n_individuals,
        haplotype_freq_rmse=_haplotype_freq_rmse(stats),
        observed_homozygosity=stats.observed_homozygosity,
        expected_homozygosity_m_based=stats.expected_homozygosity_m_based,
        expected_homozygosity_pooled_hwe=stats.expected_homozygosity_pooled_hwe,
        tv_distance_vs_m_based=stats.tv_distance_vs_m_based,
        tv_distance_vs_pooled_hwe=stats.tv_distance_vs_pooled_hwe,
        origin_pair_proportions=dict(stats.origin_pair_proportions),
        transform_enabled=transform_enabled,
        transform_reached=(tresult.reached if tresult is not None else None),
        transform_achieved_homozygosity=(
            tresult.achieved_homozygosity if tresult is not None else None
        ),
        transform_reason=(tresult.reason if tresult is not None else None),
    )


@dataclass(slots=True)
class SimulationResult:
    config: SimulationConfig
    individuals: list[SimulatedIndividual]
    statistics: SummaryStatistics
    manifest: dict = field(default_factory=dict)
    transform_results: list[TransformationResult] = field(default_factory=list)
    original_individuals: list[SimulatedIndividual] = field(default_factory=list)
    replicate_summary: ReplicateBatchSummary | None = None

    def write(self, output_dir: str | Path | None = None) -> dict[str, str]:
        out = Path(output_dir) if output_dir is not None else Path(self.config.output_dir)
        out.mkdir(parents=True, exist_ok=True)
        compress = self.config.compress

        files: dict[str, str] = {}
        files["diplotypes"] = str(
            write_diplotypes(self.individuals, out / "diplotypes.tsv", compress)
        )
        files["haplotype_frequencies"] = str(
            write_haplotype_frequencies(
                self.statistics, out / "haplotype_frequencies.tsv", compress
            )
        )
        files["diplotype_frequencies"] = str(
            write_diplotype_frequencies(
                self.statistics, out / "diplotype_frequencies.tsv", compress
            )
        )
        files["summary"] = str(write_summary(self.statistics, out / "summary.json"))

        # Machine-readable per-replicate + aggregate statistics.
        if self.replicate_summary is not None:
            rpath = out / "replicates.json"
            rpath.write_text(
                json.dumps(self.replicate_summary.to_dict(), indent=2), encoding="utf-8"
            )
            files["replicates"] = str(rpath)

        # When a transform was applied, also write the original (pre-transform)
        # diplotypes so the untransformed sample remains recoverable.
        if self.transform_results and self.original_individuals:
            files["diplotypes_original"] = str(
                write_diplotypes(
                    self.original_individuals, out / "diplotypes_original.tsv", compress
                )
            )

        manifest = dict(self.manifest)
        manifest["generated_files"] = files
        files["manifest"] = str(write_manifest(manifest, out / "manifest.json"))
        return files


class Simulator:
    """Orchestrates loading, validation, preprocessing, sampling, and statistics."""

    def __init__(self, config: SimulationConfig) -> None:
        self.config = config
        self.source_populations: list[SourcePopulation] = []
        self._renormalizations: dict[str, Renormalization] = {}
        self._aggregations: dict[str, AggregationProvenance] = {}
        self._input_paths: list[Path] = []

    def load_populations(self) -> list[SourcePopulation]:
        opts = HFXOptions(lenient=self.config.lenient_hfx)
        pops: list[SourcePopulation] = []
        self._input_paths = []
        for spec, path in zip(
            self.config.populations, self.config.resolved_hfx_paths(), strict=True
        ):
            self._input_paths.append(Path(path))
            pop = read_hfx(path, opts)
            pop = SourcePopulation(
                population_id=spec.id,
                distribution=pop.distribution,
                label=pop.label,
                provenance=pop.provenance,
            )
            pops.append(pop)
        self.source_populations = pops
        return pops

    def _preprocess_and_validate(self) -> None:
        """Apply explicit preprocessing (aggregation, authorized renorm) and validate."""
        vopts = ValidationOptions(
            allow_renormalization=self.config.allow_renormalization,
            renorm_tolerance=self.config.renorm_tolerance,
        )
        processed: list[SourcePopulation] = []
        for pop in self.source_populations:
            dist = pop.distribution

            # (rev 2) Explicit duplicate aggregation preprocessing, opt-in only.
            if self.config.aggregate_duplicates:
                dist, agg = aggregate_duplicates(dist)
                if agg.n_groups_merged:
                    self._aggregations[pop.population_id] = agg

            # Validate (strict by default). Duplicates are errors unless aggregated.
            report = validate_distribution(dist, vopts)
            report.raise_if_errors()

            # (rev 1) Authorized renormalization, recorded.
            total = dist.total_frequency()
            if self.config.allow_renormalization and abs(total - 1.0) > vopts.sum_tolerance:
                dist, renorm = renormalize(dist)
                self._renormalizations[pop.population_id] = renorm

            processed.append(
                SourcePopulation(
                    population_id=pop.population_id,
                    distribution=dist,
                    label=pop.label,
                    provenance=pop.provenance,
                )
            )
        self.source_populations = processed

    def _check_locus_consistency(self) -> None:
        if len(self.source_populations) <= 1:
            return
        ref = self.source_populations[0].distribution.locus_names
        for pop in self.source_populations[1:]:
            if pop.distribution.locus_names != ref:
                raise ValueError(
                    f"locus set mismatch: {pop.population_id!r} has "
                    f"{pop.distribution.locus_names} but expected {ref}"
                )

    def validate(self) -> None:
        """Load, preprocess, and validate all inputs; raise on error."""
        if not self.source_populations:
            self.load_populations()
        self._check_locus_consistency()
        self._preprocess_and_validate()

    def _replicate_seed_sequences(self, rep: int) -> tuple:
        """Return (sampling_ss, transformation_ss) for replicate index ``rep``.

        Uses ``SeedSequence.spawn`` whose children are deterministic by index
        and independent of how many are spawned, so replicate ``rep`` is
        identical whether generated alone or inside any larger batch.
        """
        child = np.random.SeedSequence(self.config.seed).spawn(rep + 1)[rep]
        sampling_ss, transformation_ss = child.spawn(2)
        return sampling_ss, transformation_ss

    def generate_replicate(
        self, rep: int, origin, sampler: Sampler, id_prefix: str
    ) -> tuple[list[SimulatedIndividual], list[SimulatedIndividual], object]:
        """Generate one replicate: (sampled_original, final, transform_result|None)."""
        config = self.config
        sampling_ss, transformation_ss = self._replicate_seed_sequences(rep)
        sampling_seed = int(sampling_ss.generate_state(1)[0])
        rng, _ = make_rng(sampling_seed)
        sampled = list(
            sampler.sample(config.n_individuals, rng, id_prefix=id_prefix, replicate=rep)
        )
        tcfg = config.transformation
        if tcfg.enabled:
            t_seed = int(transformation_ss.generate_state(1)[0])
            t_rng, _ = make_rng(t_seed)
            target = TransformationTarget(
                fraction=tcfg.target_fraction, delta=tcfg.target_delta
            )
            final, tresult = transform_homozygosity(sampled, target, rng=t_rng)
            return sampled, final, tresult
        return sampled, sampled, None

    def run(self) -> SimulationResult:
        self.validate()
        config = self.config
        origin = _build_origin_model(config, len(self.source_populations))
        sampler = Sampler(self.source_populations, origin)

        canonical = _canonical_config(config)
        software_version = _software_version()
        run_hash = _run_hash(canonical, software_version)
        id_prefix = f"DF-{run_hash}"

        tcfg = config.transformation
        original_individuals: list[SimulatedIndividual] = []
        all_individuals: list[SimulatedIndividual] = []
        transform_results = []
        replicate_records: list[ReplicateStatistics] = []
        generator_name = "PCG64"
        for rep in range(config.replicates):
            sampled, rep_final, tresult = self.generate_replicate(
                rep, origin, sampler, id_prefix
            )
            original_individuals.extend(sampled)
            if tresult is not None:
                transform_results.append(tresult)
            all_individuals.extend(rep_final)

            # Per-replicate statistics on the replicate's FINAL sample.
            rep_stats = compute_summary(rep_final, self.source_populations, origin.matrix)
            replicate_records.append(
                _replicate_record(rep, rep_stats, tcfg.enabled, tresult)
            )

        batch_summary = aggregate_replicates(replicate_records)

        # Statistics on the FINAL (possibly transformed) sample.
        stats = compute_summary(all_individuals, self.source_populations, origin.matrix)
        self._transform_results = transform_results
        self._original_individuals = original_individuals

        manifest = {
            "software_version": software_version,
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "created": datetime.now(UTC).isoformat(),
            "run_hash": run_hash,
            "canonical_config": canonical,
            "rng": {
                "generator": generator_name,
                "seed": config.seed,
                "replicate_seed_derivation": "numpy.random.SeedSequence.spawn",
                "stream_hierarchy": (
                    "root.spawn(replicates); each replicate .spawn(2) -> "
                    "[sampling, transformation]. Sampling is derived independently "
                    "of transformation, so enabling the transform does not change "
                    "the sampled diplotypes."
                ),
            },
            "reproducibility": (
                "Byte-identical output is guaranteed for identical inputs, seed, "
                "DiplotypeForge version, canonical configuration, AND numpy version. "
                "Output is NOT guaranteed identical across different numpy versions; "
                "pin the recorded numpy_version to reproduce exactly."
            ),
            "inputs": [
                {
                    "population_id": pop.population_id,
                    "source_path": str(path),
                    "sha256": _sha256(path) if path.exists() else None,
                }
                for pop, path in zip(self.source_populations, self._input_paths, strict=False)
            ],
            "hfx_schema": schema_provenance(),
            "populations": [
                {
                    "id": p.population_id,
                    "label": p.label,
                    "n_haplotypes": p.distribution.n_entries,
                    "input_frequency_sum": p.distribution.total_frequency(),
                }
                for p in self.source_populations
            ],
            "origin_model": {
                "mode": config.mode,
                "matrix": origin.matrix.tolist(),
                "symmetric": bool(np.allclose(origin.matrix, origin.matrix.T)),
                "population_order": [p.population_id for p in self.source_populations],
                "marginal_ancestry_proportions": dict(
                    zip(
                        [p.population_id for p in self.source_populations],
                        marginal_ancestry_proportions(origin.matrix).tolist(),
                        strict=True,
                    )
                ),
            },
        }
        if self._renormalizations:
            manifest["renormalization"] = {
                pid: {
                    "original_sum": r.original_sum,
                    "scale_factor": r.scale_factor,
                    "deviation": r.tolerance,
                }
                for pid, r in self._renormalizations.items()
            }
        if self._aggregations:
            manifest["duplicate_aggregation"] = {
                pid: agg.merged for pid, agg in self._aggregations.items()
            }

        if tcfg.enabled and transform_results:
            manifest["transformation"] = {
                "enabled": True,
                "type": tcfg.type,
                "strategy": tcfg.strategy,
                "target_fraction": tcfg.target_fraction,
                "target_delta": tcfg.target_delta,
                "per_replicate": [asdict(tr) for tr in transform_results],
            }
            # Diagnostics on the first replicate (before vs after vs M vs pooled).
            n_per = config.n_individuals
            before0 = original_individuals[:n_per]
            after0 = all_individuals[:n_per]
            diag = transform_diagnostics(
                before0, after0, self.source_populations, origin.matrix
            )
            manifest["transformation"]["diagnostics_replicate_0"] = {
                "observed_homozygosity_before": diag.observed_homozygosity_before,
                "observed_homozygosity_after": diag.observed_homozygosity_after,
                "tv_before_vs_m": diag.tv_before_vs_m,
                "tv_after_vs_m": diag.tv_after_vs_m,
                "tv_before_vs_pooled_hwe": diag.tv_before_vs_pooled,
                "tv_after_vs_pooled_hwe": diag.tv_after_vs_pooled,
                "tv_before_vs_after": diag.tv_before_vs_after,
                "origin_pair_before": diag.origin_pair_before,
                "origin_pair_after": diag.origin_pair_after,
                "note": diag.note,
            }

        manifest["replicates"] = {
            "n_replicates": batch_summary.n_replicates,
            "n_success": batch_summary.n_success,
            "n_fail": batch_summary.n_fail,
            "aggregates": batch_summary.aggregates,
        }

        return SimulationResult(
            config=config,
            individuals=all_individuals,
            statistics=stats,
            manifest=manifest,
            transform_results=transform_results,
            original_individuals=original_individuals,
            replicate_summary=batch_summary,
        )

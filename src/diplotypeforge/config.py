"""Simulation configuration (TOML-backed).

Milestone 1 supports the ``single_population`` mode. The config schema is
forward-looking: it accepts ``mode`` and optional ``matrix``/``weights`` blocks
so Milestone 2 modes slot in without breaking changes.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

SUPPORTED_MODES = frozenset(
    {"single_population", "wahlund", "admixture", "matrix"}
)
# Modes implemented and enabled (Milestone 2 enables all four).
ENABLED_MODES = frozenset({"single_population", "wahlund", "admixture", "matrix"})


class ConfigError(ValueError):
    """Raised for malformed or unsupported simulation configuration."""


@dataclass(frozen=True, slots=True)
class PopulationSpec:
    id: str
    hfx: str


TRANSFORM_TYPES = frozenset({"homozygosity"})
TRANSFORM_STRATEGIES = frozenset({"reassignment"})


@dataclass(frozen=True, slots=True)
class TransformationConfig:
    """The ``[transformation]`` block (Milestone 3)."""

    enabled: bool = False
    type: str = "homozygosity"
    strategy: str = "reassignment"
    target_fraction: float | None = None
    target_delta: float | None = None

    def __post_init__(self) -> None:
        if not self.enabled:
            return
        if self.type not in TRANSFORM_TYPES:
            raise ConfigError(
                f"unknown transformation type {self.type!r}; supported: {sorted(TRANSFORM_TYPES)}"
            )
        if self.strategy not in TRANSFORM_STRATEGIES:
            raise ConfigError(
                f"unknown transformation strategy {self.strategy!r}; "
                f"supported: {sorted(TRANSFORM_STRATEGIES)}"
            )
        if (self.target_fraction is None) == (self.target_delta is None):
            raise ConfigError(
                "transformation requires exactly one of target_fraction or target_delta"
            )
        if self.target_fraction is not None and not 0.0 <= self.target_fraction <= 1.0:
            raise ConfigError(
                f"target_fraction must be in [0, 1]; got {self.target_fraction}"
            )


@dataclass(frozen=True, slots=True)
class SimulationConfig:
    n_individuals: int
    seed: int
    mode: str = "single_population"
    output_dir: str = "output"
    compress: bool = False
    replicates: int = 1
    strict: bool = True
    lenient_hfx: bool = False
    allow_renormalization: bool = False
    renorm_tolerance: float = 0.0
    aggregate_duplicates: bool = False
    populations: tuple[PopulationSpec, ...] = field(default_factory=tuple)
    # Forward-looking mode parameters (Milestone 2):
    admixture_a: float | None = None
    weights: tuple[float, ...] | None = None
    matrix: tuple[tuple[float, ...], ...] | None = None
    symmetric: bool = False
    # Milestone 3:
    transformation: TransformationConfig = field(default_factory=TransformationConfig)
    base_dir: str = "."

    def __post_init__(self) -> None:
        if self.n_individuals < 0:
            raise ConfigError("n_individuals must be non-negative")
        if self.mode not in SUPPORTED_MODES:
            raise ConfigError(
                f"unknown mode {self.mode!r}; supported: {sorted(SUPPORTED_MODES)}"
            )
        if self.mode not in ENABLED_MODES:
            raise ConfigError(
                f"mode {self.mode!r} is not enabled (available: {sorted(ENABLED_MODES)})"
            )
        if not self.populations:
            raise ConfigError("at least one [[populations]] entry is required")
        ids = [p.id for p in self.populations]
        if len(ids) != len(set(ids)):
            raise ConfigError(f"duplicate population ids: {ids}")
        if self.replicates < 1:
            raise ConfigError("replicates must be >= 1")

        n_pops = len(self.populations)
        if self.mode == "single_population":
            if n_pops != 1:
                raise ConfigError(
                    f"single_population mode requires exactly one population, got {n_pops}"
                )
        elif self.mode == "wahlund":
            if n_pops < 2:
                raise ConfigError(f"wahlund mode requires >= 2 populations, got {n_pops}")
            if self.weights is not None and len(self.weights) != n_pops:
                raise ConfigError(
                    f"wahlund weights length {len(self.weights)} != {n_pops} populations"
                )
        elif self.mode == "admixture":
            if n_pops != 2:
                raise ConfigError(
                    f"admixture mode is the balanced two-population model; requires exactly "
                    f"2 populations, got {n_pops}"
                )
            if self.admixture_a is None:
                raise ConfigError("admixture mode requires 'admixture_a' in [0, 1]")
            if not 0.0 <= float(self.admixture_a) <= 1.0:
                raise ConfigError(f"admixture_a must be in [0, 1]; got {self.admixture_a}")
        elif self.mode == "matrix":
            if self.matrix is None:
                raise ConfigError("matrix mode requires a 'matrix' (K x K) in [simulation]")
            k = len(self.matrix)
            if k != n_pops:
                raise ConfigError(
                    f"matrix is {k}x? but there are {n_pops} populations; dims must match"
                )
            for row in self.matrix:
                if len(row) != k:
                    raise ConfigError(f"matrix must be square KxK; got a row of length {len(row)}")

    def resolved_hfx_paths(self) -> list[Path]:
        base = Path(self.base_dir)
        return [(base / p.hfx) if not Path(p.hfx).is_absolute() else Path(p.hfx)
                for p in self.populations]

    @classmethod
    def from_dict(cls, data: dict, base_dir: str = ".") -> SimulationConfig:
        sim = data.get("simulation", {})
        pops_raw = data.get("populations", [])
        if not pops_raw:
            raise ConfigError("configuration must define [[populations]]")
        populations = tuple(
            PopulationSpec(id=str(p["id"]), hfx=str(p["hfx"])) for p in pops_raw
        )
        weights = sim.get("weights")
        matrix = sim.get("matrix")
        tr = data.get("transformation", {})
        transformation = TransformationConfig(
            enabled=bool(tr.get("enabled", False)),
            type=str(tr.get("type", "homozygosity")),
            strategy=str(tr.get("strategy", "reassignment")),
            target_fraction=tr.get("target_fraction"),
            target_delta=tr.get("target_delta"),
        )
        return cls(
            n_individuals=int(sim.get("n_individuals", 0)),
            seed=int(sim.get("seed", 0)),
            mode=str(sim.get("mode", "single_population")),
            output_dir=str(sim.get("output_dir", "output")),
            compress=bool(sim.get("compress", False)),
            replicates=int(sim.get("replicates", 1)),
            strict=bool(sim.get("strict", True)),
            lenient_hfx=bool(sim.get("lenient_hfx", False)),
            allow_renormalization=bool(sim.get("allow_renormalization", False)),
            renorm_tolerance=float(sim.get("renorm_tolerance", 0.0)),
            aggregate_duplicates=bool(sim.get("aggregate_duplicates", False)),
            populations=populations,
            admixture_a=sim.get("admixture_a"),
            weights=tuple(weights) if weights is not None else None,
            matrix=tuple(tuple(row) for row in matrix) if matrix is not None else None,
            symmetric=bool(sim.get("symmetric", False)),
            transformation=transformation,
            base_dir=base_dir,
        )

    @classmethod
    def from_file(cls, path: str | Path) -> SimulationConfig:
        p = Path(path)
        if not p.exists():
            raise ConfigError(f"configuration file not found: {p}")
        try:
            data = tomllib.loads(p.read_text(encoding="utf-8"))
        except tomllib.TOMLDecodeError as exc:
            raise ConfigError(f"invalid TOML in {p}: {exc}") from exc
        return cls.from_dict(data, base_dir=str(p.parent))

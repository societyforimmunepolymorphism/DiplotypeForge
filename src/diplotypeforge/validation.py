"""Validation and preprocessing of haplotype-frequency distributions.

Default posture is **strict** (reviewer revisions 1 & 2):

  * Frequency sum must equal 1.0 within a tiny round-off ``sum_tolerance``.
    A larger renormalization is applied only when explicitly authorized
    (``allow_renormalization`` + ``renorm_tolerance``) and is recorded.
  * Duplicate haplotype identifiers are **errors**. Aggregating duplicates is an
    explicit preprocessing operation (:func:`aggregate_duplicates`) that warns
    and records provenance; it is never implicit.
  * Negative / missing frequency → error. Zero frequency → warning (unsamplable).

Returns a structured :class:`ValidationReport`; ``raise_if_errors`` is provided.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, replace
from enum import StrEnum

from .model.population import FrequencyEntry, HaplotypeFrequencyDistribution


class Severity(StrEnum):
    ERROR = "error"
    WARNING = "warning"


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    severity: Severity
    code: str
    message: str


@dataclass(slots=True)
class ValidationReport:
    issues: list[ValidationIssue] = field(default_factory=list)

    def add(self, severity: Severity, code: str, message: str) -> None:
        self.issues.append(ValidationIssue(severity, code, message))

    @property
    def errors(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity is Severity.ERROR]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity is Severity.WARNING]

    @property
    def ok(self) -> bool:
        return not self.errors

    def raise_if_errors(self) -> None:
        if self.errors:
            lines = "\n".join(f"  [{i.code}] {i.message}" for i in self.errors)
            raise ValidationError(f"frequency distribution validation failed:\n{lines}")


class ValidationError(ValueError):
    """Raised by :meth:`ValidationReport.raise_if_errors` when errors exist."""


@dataclass(frozen=True, slots=True)
class ValidationOptions:
    """Tunable validation behavior. Defaults are strict (see OPEN_QUESTIONS Q1/Q2).

    ``sum_tolerance`` absorbs floating-point round-off only. Any deviation beyond
    it is an error unless ``allow_renormalization`` is set, in which case a
    deviation up to ``renorm_tolerance`` is permitted (and recorded by the
    caller via :func:`renormalize`).
    """

    sum_tolerance: float = 1e-6
    allow_renormalization: bool = False
    renorm_tolerance: float = 0.0
    allow_zero_frequency: bool = True


def validate_distribution(
    dist: HaplotypeFrequencyDistribution,
    options: ValidationOptions | None = None,
) -> ValidationReport:
    """Validate a distribution, returning a structured report."""
    opts = options or ValidationOptions()
    report = ValidationReport()

    n_loci = len(dist.loci)
    locus_names = dist.locus_names

    if dist.n_entries == 0:
        report.add(Severity.ERROR, "empty", "distribution has no haplotype entries")
        return report

    seen: dict[tuple[str, ...], int] = {}
    for idx, entry in enumerate(dist.entries):
        hap = entry.haplotype
        freq = entry.frequency

        if hap.n_loci != n_loci:
            report.add(
                Severity.ERROR,
                "inconsistent-loci",
                f"entry {idx}: haplotype spans {hap.n_loci} loci but distribution "
                f"declares {n_loci}",
            )
        else:
            for allele, expected in zip(hap.alleles, locus_names, strict=True):
                if allele.locus != expected:
                    report.add(
                        Severity.ERROR,
                        "invalid-allele",
                        f"entry {idx}: allele {allele.name!r} has locus "
                        f"{allele.locus!r} but expected {expected!r}",
                    )

        if freq is None:  # type: ignore[comparison-overlap]
            report.add(Severity.ERROR, "missing-frequency", f"entry {idx}: missing frequency")
        elif freq < 0:
            report.add(
                Severity.ERROR, "negative-frequency", f"entry {idx}: negative frequency {freq}"
            )
        elif freq == 0:
            sev = Severity.WARNING if opts.allow_zero_frequency else Severity.ERROR
            report.add(sev, "zero-frequency", f"entry {idx}: zero frequency (cannot be sampled)")

        # Duplicate haplotype identifiers are ERRORS by default (rev 2).
        key = hap.key()
        if key in seen:
            report.add(
                Severity.ERROR,
                "duplicate-haplotype",
                f"entry {idx}: duplicate haplotype also at entry {seen[key]}: "
                f"{'~'.join(key)}. Use aggregate_duplicates preprocessing to merge.",
            )
        else:
            seen[key] = idx

    # --- frequency sum ------------------------------------------------------
    total = dist.total_frequency()
    deviation = abs(total - 1.0)
    if deviation <= opts.sum_tolerance:
        pass  # pure round-off; acceptable
    elif opts.allow_renormalization and deviation <= opts.renorm_tolerance:
        report.add(
            Severity.WARNING,
            "sum-renormalized",
            f"frequency sum {total:.6g} deviates by {deviation:.3g}; authorized "
            f"renormalization to 1.0 will be applied and recorded",
        )
    else:
        hint = (
            "set allow_renormalization=true with an adequate renorm_tolerance to permit it"
            if not opts.allow_renormalization
            else f"deviation exceeds renorm_tolerance {opts.renorm_tolerance:g}"
        )
        report.add(
            Severity.ERROR,
            "sum-out-of-tolerance",
            f"frequency sum {total:.6g} deviates from 1.0 by {deviation:.3g} "
            f"(round-off tolerance {opts.sum_tolerance:g}); {hint}",
        )

    return report


# --- explicit preprocessing operations --------------------------------------


@dataclass(frozen=True, slots=True)
class Renormalization:
    """Provenance record of an applied renormalization."""

    original_sum: float
    scale_factor: float
    tolerance: float


def renormalize(
    dist: HaplotypeFrequencyDistribution,
) -> tuple[HaplotypeFrequencyDistribution, Renormalization]:
    """Scale all frequencies so they sum to 1.0. Returns the new distribution
    and a provenance record. The caller is responsible for authorization.
    """
    total = dist.total_frequency()
    if total <= 0:
        raise ValidationError("cannot renormalize a distribution with non-positive total")
    scale = 1.0 / total
    new_entries = tuple(
        replace(e, frequency=e.frequency * scale) for e in dist.entries
    )
    new_dist = replace(dist, entries=new_entries)
    return new_dist, Renormalization(
        original_sum=total, scale_factor=scale, tolerance=abs(total - 1.0)
    )


@dataclass(frozen=True, slots=True)
class AggregationProvenance:
    """Provenance of an explicit duplicate-aggregation preprocessing step."""

    merged: dict[str, list[float]] = field(default_factory=dict)

    @property
    def n_groups_merged(self) -> int:
        return len(self.merged)


def aggregate_duplicates(
    dist: HaplotypeFrequencyDistribution,
) -> tuple[HaplotypeFrequencyDistribution, AggregationProvenance]:
    """Explicitly merge duplicate haplotypes by summing their frequencies.

    Emits a warning and records provenance (which keys merged and their
    contributing frequencies). This is an opt-in preprocessing step; it is never
    applied implicitly by validation or the adapter.
    """
    order: list[tuple[str, ...]] = []
    acc: dict[tuple[str, ...], FrequencyEntry] = {}
    contributions: dict[tuple[str, ...], list[float]] = {}

    for entry in dist.entries:
        key = entry.haplotype.key()
        contributions.setdefault(key, []).append(entry.frequency)
        if key in acc:
            prev = acc[key]
            merged_count = None
            if prev.count is not None or entry.count is not None:
                merged_count = (prev.count or 0.0) + (entry.count or 0.0)
            acc[key] = replace(prev, frequency=prev.frequency + entry.frequency, count=merged_count)
        else:
            acc[key] = entry
            order.append(key)

    merged = {
        "~".join(k): v for k, v in contributions.items() if len(v) > 1
    }
    if merged:
        warnings.warn(
            f"aggregate_duplicates merged {len(merged)} duplicate haplotype group(s): "
            f"{list(merged)}",
            stacklevel=2,
        )
    new_dist = replace(dist, entries=tuple(acc[k] for k in order))
    return new_dist, AggregationProvenance(merged=merged)

"""Fit stutter tables from aggregate homopolymer targets and validate simulations.

Used by ``helpers/calibrate_read_profile.py`` (empirical engine). Inputs are the
``hp_P_obs_given_true`` aggregates: per key (``C8|+``) the number of observed
runs ``n`` and ``p_obs_minus_true``, the pmf of observed minus true run length.

With the empirical error channel homopolymer runs receive no base-level errors,
so a fitted entry equals the measured pmf, truncated to ``|delta| <= max_delta``
(and to deltas that keep at least zero bases) and renormalised. The fallback
for unfitted runs (see :mod:`.stutter`) is derived from the same fitted table:

* ``log_odds_slope_per_base``: for every base/strand with at least
  ``min_lengths_for_slope`` fitted lengths, the least-squares slope of the
  error log-odds against run length; the median over those groups, clamped at
  0 so extrapolation never makes longer runs more accurate;
* ``max_extrapolation_bases``: configured cap on how far the slope is
  extrapolated beyond the fitted range;
* ``generic``: the ``n``-weighted pool of the fitted entries at
  ``generic_ref_len`` over all bases and strands.

Pure functions, no I/O.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from collections.abc import Mapping
from typing import Any

from .stutter import StutterFallback

ROUND_DIGITS = 5


def _key_parts(key: str) -> tuple[str, int, str]:
    run, _, strand = key.partition("|")
    return run[0], int(run[1:]), strand


def _target_pmf(entry: Mapping[str, Any]) -> dict[int, float]:
    return {int(d): float(p) for d, p in entry["p_obs_minus_true"].items()}


def _normalise(pmf: Mapping[int, float], name: str) -> dict[int, float]:
    total = sum(p for p in pmf.values() if p > 0)
    if total <= 0:
        raise ValueError(f"stutter pmf for '{name}' has no probability mass to normalise")
    return {d: p / total for d, p in sorted(pmf.items()) if p > 0}


def fit_stutter_pmfs(
    hp_targets: Mapping[str, Mapping[str, Any]], *, min_n: int, min_len: int, max_delta: int
) -> dict[str, dict[int, float]]:
    """Per-strand stutter pmfs for every target key with ``n >= min_n`` and length >= ``min_len``."""
    pmfs: dict[str, dict[int, float]] = {}
    for key, entry in sorted(hp_targets.items()):
        _, length, strand = _key_parts(key)
        if strand == "both" or length < min_len or entry["n"] < min_n:
            continue
        lower = -min(length, max_delta)
        kept = {d: p for d, p in _target_pmf(entry).items() if lower <= d <= max_delta}
        pmfs[key] = _normalise(kept, f"{key} within |delta| <= {max_delta}")
    return pmfs


def _log_odds(pmf: Mapping[int, float]) -> float | None:
    correct = pmf.get(0, 0.0)
    if correct <= 0 or correct >= 1:
        return None
    return math.log((1.0 - correct) / correct)


def fit_log_odds_slope(pmfs: Mapping[str, Mapping[int, float]], *, min_lengths: int) -> float:
    """Median least-squares slope of error log-odds vs run length, clamped at 0."""
    groups: dict[tuple[str, str], list[tuple[int, float]]] = defaultdict(list)
    for key, pmf in pmfs.items():
        base, length, strand = _key_parts(key)
        log_odds = _log_odds(pmf)
        if log_odds is not None:
            groups[(base, strand)].append((length, log_odds))
    slopes = []
    for points in groups.values():
        if len(points) < min_lengths:
            continue
        mean_x = statistics.fmean(x for x, _ in points)
        mean_y = statistics.fmean(y for _, y in points)
        sxx = sum((x - mean_x) ** 2 for x, _ in points)
        sxy = sum((x - mean_x) * (y - mean_y) for x, y in points)
        slopes.append(sxy / sxx)
    if not slopes:
        raise ValueError(
            f"no base/strand has >= min_lengths ({min_lengths}) fitted run lengths; "
            "lower min_lengths or set the slope explicitly"
        )
    return max(0.0, statistics.median(slopes))


def pooled_generic_pmf(
    hp_targets: Mapping[str, Mapping[str, Any]],
    pmfs: Mapping[str, Mapping[int, float]],
    *,
    ref_len: int,
) -> dict[int, float]:
    """``n``-weighted pool of the fitted pmfs at ``ref_len`` (all bases and strands)."""
    pooled: dict[int, float] = defaultdict(float)
    for key, pmf in pmfs.items():
        if _key_parts(key)[1] != ref_len:
            continue
        for delta, p in pmf.items():
            pooled[delta] += hp_targets[key]["n"] * p
    if not pooled:
        raise ValueError(f"no fitted stutter entry at generic ref_len {ref_len}")
    return _normalise(pooled, f"generic pool at length {ref_len}")


def round_pmf(pmf: Mapping[int, float]) -> dict[str, float]:
    """Round to ``ROUND_DIGITS`` so the pmf still sums to 1.

    The rounding drift goes to delta 0, or to the most probable delta when
    that would make P(0) negative.
    """
    rounded = {d: round(p, ROUND_DIGITS) for d, p in sorted(pmf.items())}
    rounded = {d: p for d, p in rounded.items() if p > 0}
    drift = 1.0 - sum(rounded.values())
    target = 0 if rounded.get(0, 0.0) + drift >= 0 else max(rounded, key=lambda d: rounded[d])
    rounded[target] = round(rounded.get(target, 0.0) + drift, ROUND_DIGITS)
    return {str(d): p for d, p in sorted(rounded.items()) if p > 0}


def stutter_profile_section(
    hp_targets: Mapping[str, Mapping[str, Any]],
    *,
    min_n: int,
    min_len: int,
    max_delta: int,
    min_lengths_for_slope: int,
    generic_ref_len: int,
    max_extrapolation_bases: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return (molecule-model stutter keys, provenance record) for a read profile."""
    pmfs = fit_stutter_pmfs(hp_targets, min_n=min_n, min_len=min_len, max_delta=max_delta)
    slope = fit_log_odds_slope(pmfs, min_lengths=min_lengths_for_slope)
    generic = pooled_generic_pmf(hp_targets, pmfs, ref_len=generic_ref_len)
    fallback = {
        "rule": "log_odds_interpolate",
        "log_odds_slope_per_base": round(slope, ROUND_DIGITS),
        "max_extrapolation_bases": max_extrapolation_bases,
        "generic": {"ref_len": generic_ref_len, "pmf": round_pmf(generic)},
    }
    StutterFallback.from_dict(fallback)  # validate before writing
    per_strand = sorted(k for k in hp_targets if _key_parts(k)[2] != "both")
    provenance = {
        "method": (
            "empirical engine: fitted entry = measured hp_P_obs_given_true pmf truncated to "
            "|delta| <= max_delta (and >= -run length) and renormalised"
        ),
        "min_target_n": min_n,
        "max_extrapolation_bases": (
            f"{max_extrapolation_bases}: the log-odds slope is only extrapolated this many "
            "bases beyond the fitted range; longer runs are treated like the capped length"
        ),
        "min_run_len": min_len,
        "max_delta": max_delta,
        "fitted_keys": sorted(pmfs),
        "excluded_keys": {
            key: f"n={hp_targets[key]['n']}"
            for key in per_strand
            if key not in pmfs and _key_parts(key)[1] >= min_len
        },
        "rule": (
            "log_odds_interpolate: runs without a fitted entry of the same base and strand "
            "(strand 'both' if no strand-specific entry) interpolate the error log-odds and "
            "error shape linearly between the nearest fitted lengths on both sides; beyond "
            "the longest (or below the shortest) fitted length the error odds of that entry "
            "are scaled by exp(log_odds_slope_per_base * length difference), with the "
            "difference capped at max_extrapolation_bases; bases without "
            "fitted entries, or references without a usable error delta, use the generic pmf"
        ),
        "slope": (
            f"median least-squares slope of error log-odds vs length over base/strand "
            f"groups with >= {min_lengths_for_slope} fitted lengths, clamped at 0"
        ),
        "generic": (
            f"n-weighted pool of fitted entries at length {generic_ref_len}; used with a "
            "warning for bases without any fitted entry"
        ),
    }
    stutter = {key: round_pmf(pmf) for key, pmf in pmfs.items()}
    return {"stutter": stutter, "stutter_fallback": fallback}, provenance


def compare_to_targets(
    observed: Mapping[str, Counter[int]],
    hp_targets: Mapping[str, Mapping[str, Any]],
    *,
    max_delta: int,
    tolerance: float,
) -> dict[str, dict[str, Any]]:
    """Per-strand observed-vs-target homopolymer spectra for every target key.

    ``max_abs_residual`` is the largest absolute pmf difference over deltas in
    ``[-max_delta, max_delta]``; a key passes when it is within ``tolerance``.
    Keys whose run does not occur in the simulated template are reported as
    ``not_in_template`` with ``within_tolerance`` None.
    """
    report: dict[str, dict[str, Any]] = {}
    for key, entry in sorted(hp_targets.items()):
        if _key_parts(key)[2] == "both":
            continue
        target = _target_pmf(entry)
        counts = observed.get(key, Counter())
        n_sim = sum(counts.values())
        row: dict[str, Any] = {
            "n_target": entry["n"],
            "n_sim": n_sim,
            "p_correct_target": round(target.get(0, 0.0), 4),
        }
        if n_sim == 0:
            row.update(status="not_in_template", within_tolerance=None)
        else:
            sim = {d: c / n_sim for d, c in counts.items()}
            residual = max(
                abs(sim.get(d, 0.0) - target.get(d, 0.0)) for d in range(-max_delta, max_delta + 1)
            )
            row.update(
                status="compared",
                p_correct_sim=round(sim.get(0, 0.0), 4),
                max_abs_residual=round(residual, 4),
                within_tolerance=residual <= tolerance,
            )
        report[key] = row
    return report

"""Strand-aware homopolymer stutter tables and their fallback for unfitted runs.

A :class:`StutterTable` maps keys such as ``"C7|+"`` (base and run length in
source orientation, read strand ``+``/``-`` or ``both``) to a probability mass
over length deltas. Fitted entries come from measured data.

Runs without a fitted entry are resolved by an optional
:class:`StutterFallback` (rule ``log_odds_interpolate``). Error log-odds are
``log(P(delta != 0) / P(delta == 0))``; the error shape is the pmf conditional
on an error. Deltas that would remove more bases than the run has are dropped.

1. **Fitted lengths.** Use the fitted entries of the same base and strand (or
   of strand ``both`` when that base has no strand-specific entries).
2. **Interpolation.** A run between two fitted lengths takes the log-odds
   interpolated linearly between the nearest fitted lengths on both sides, and
   the error shape mixed with the same weights.
3. **Extrapolation.** A run longer than the longest (or shorter than the
   shortest) fitted length scales that entry's error odds by
   ``exp(log_odds_slope_per_base * d)`` and keeps its error shape, where ``d``
   is the length difference clamped to ``±max_extrapolation_bases``: runs
   further from the data are treated like the capped length.
4. **Generic.** A base without fitted entries uses the profile's generic pmf
   at its ``ref_len``, scaled as in 3, and logs a warning (once per base and
   table). The generic pmf is also used when the reference keeps no error
   delta for this run, so a run resolved by the fallback is never error-free;
   it must itself keep error mass at the table's minimum run length.

Without a fallback a missing key draws no delta. That is only safe when the
sequencer adds its own homopolymer errors (pbsim3); read profiles with an
empirical error channel must define a fallback (see
:func:`muc_one_up.read_simulator.read_profiles.load_read_profile`).
"""

from __future__ import annotations

import logging
import math
import random
from collections.abc import Mapping
from dataclasses import dataclass, field

Pmf = tuple[tuple[int, float], ...]

STRANDS = ("+", "-")
FALLBACK_RULES = ("log_odds_interpolate",)
_FALLBACK_KEYS = {"rule", "log_odds_slope_per_base", "max_extrapolation_bases", "generic"}
_PMF_TOLERANCE = 1e-6


def validate_stutter_key(key: str) -> None:
    """Raise ValueError unless ``key`` looks like ``C7|+``, ``G4|-`` or ``A5|both``."""
    run, sep, strand = key.partition("|")
    if (
        not sep
        or strand not in (*STRANDS, "both")
        or len(run) < 2
        or run[0] not in "ACGT"
        or not run[1:].isdigit()
    ):
        raise ValueError(f"invalid stutter key '{key}': expected e.g. 'C7|+', 'G4|-' or 'A5|both'")


def parse_pmf(data: Mapping[str, float] | Mapping[int, float], name: str, run_len: int) -> Pmf:
    """Validate a delta -> probability mapping for a run of ``run_len`` bases."""
    items = tuple(sorted((int(delta), float(p)) for delta, p in data.items()))
    if any(delta < -run_len and p > 0 for delta, p in items):
        raise ValueError(f"stutter delta for '{name}' cannot remove more than {run_len} bases")
    if any(p < 0 for _, p in items) or not math.isclose(
        sum(p for _, p in items), 1.0, abs_tol=_PMF_TOLERANCE
    ):
        raise ValueError(f"stutter pmf for '{name}' must be non-negative and sum to 1")
    return items


def scale_error_odds(pmf: Pmf, log_odds_shift: float, run_len: int) -> Pmf:
    """Shift the error log-odds of ``pmf`` and keep its error shape.

    Deltas below ``-run_len`` are dropped before scaling. A pmf without error
    mass stays error-free; one without mass at 0 stays error-only.
    """
    p_correct = dict(pmf).get(0, 0.0)
    errors = [(d, p) for d, p in pmf if d != 0 and d >= -run_len and p > 0]
    error_mass = sum(p for _, p in errors)
    if error_mass == 0:
        return ((0, 1.0),)
    if p_correct == 0:
        new_correct = 0.0
    else:
        odds = error_mass / p_correct * math.exp(log_odds_shift)
        new_correct = 1.0 / (1.0 + odds)
    scaled = [(d, (1.0 - new_correct) * p / error_mass) for d, p in errors]
    if new_correct > 0:
        scaled.append((0, new_correct))
    return tuple(sorted(scaled))


def _has_error_mass(pmf: Pmf) -> bool:
    return any(d != 0 and p > 0 for d, p in pmf)


def interpolate_log_odds(
    lower: Pmf, lower_len: int, upper: Pmf, upper_len: int, run_len: int
) -> Pmf | None:
    """Interpolate error log-odds and error shape between two fitted lengths.

    Returns None when either side has no error mass or no mass at 0 for this
    run (log-odds undefined); callers then extrapolate instead.
    """
    weight = (run_len - lower_len) / (upper_len - lower_len)
    log_odds = 0.0
    shape: dict[int, float] = {}
    for pmf, w in ((lower, 1.0 - weight), (upper, weight)):
        p_correct = dict(pmf).get(0, 0.0)
        errors = [(d, p) for d, p in pmf if d != 0 and d >= -run_len and p > 0]
        error_mass = sum(p for _, p in errors)
        if error_mass == 0 or p_correct == 0:
            return None
        log_odds += w * math.log(error_mass / p_correct)
        for d, p in errors:
            shape[d] = shape.get(d, 0.0) + w * p / error_mass
    new_correct = 1.0 / (1.0 + math.exp(log_odds))
    scaled = [(d, (1.0 - new_correct) * p) for d, p in shape.items() if p > 0]
    return tuple(sorted([*scaled, (0, new_correct)]))


@dataclass(frozen=True)
class StutterFallback:
    """Documented rule for runs without a fitted stutter entry (see module docstring)."""

    log_odds_slope_per_base: float
    max_extrapolation_bases: int
    generic_pmf: Pmf
    generic_ref_len: int
    rule: str = "log_odds_interpolate"

    def __post_init__(self) -> None:
        if self.rule not in FALLBACK_RULES:
            raise ValueError(f"stutter_fallback rule must be one of {FALLBACK_RULES}")
        slope = self.log_odds_slope_per_base
        if not math.isfinite(slope) or slope < 0:
            raise ValueError(f"log_odds_slope_per_base must be finite and >= 0, got {slope}")
        if self.max_extrapolation_bases < 0:
            raise ValueError(
                f"max_extrapolation_bases must be >= 0, got {self.max_extrapolation_bases}"
            )
        if self.generic_ref_len < 1:
            raise ValueError("stutter_fallback generic ref_len must be >= 1")
        if dict(self.generic_pmf).get(0, 0.0) >= 1.0:
            raise ValueError("stutter_fallback generic pmf must have error mass (P(0) < 1)")

    def shift(self, length_difference: int) -> float:
        """Log-odds shift for extrapolating over ``length_difference`` bases (capped)."""
        cap = self.max_extrapolation_bases
        return self.log_odds_slope_per_base * max(-cap, min(cap, length_difference))

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> StutterFallback:
        """Build from the profile's ``molecules.stutter_fallback`` mapping."""
        unknown = set(data) - _FALLBACK_KEYS
        if unknown:
            raise ValueError(f"unknown stutter_fallback keys: {sorted(unknown)}")
        generic = data["generic"]
        if not isinstance(generic, Mapping) or set(generic) != {"ref_len", "pmf"}:
            raise ValueError("stutter_fallback generic must be {'ref_len': int, 'pmf': {...}}")
        ref_len = int(generic["ref_len"])
        if ref_len < 1:
            raise ValueError("stutter_fallback generic ref_len must be >= 1")
        return cls(
            log_odds_slope_per_base=float(data["log_odds_slope_per_base"]),  # type: ignore[arg-type]
            max_extrapolation_bases=int(data["max_extrapolation_bases"]),  # type: ignore[call-overload]
            generic_pmf=parse_pmf(generic["pmf"], "stutter_fallback generic", ref_len),
            generic_ref_len=ref_len,
            rule=str(data["rule"]),
        )


@dataclass(frozen=True)
class StutterTable:
    """Per-(base, length, strand) probability mass over homopolymer length deltas."""

    pmfs: Mapping[str, Pmf]
    min_len: int = 3
    fallback: StutterFallback | None = None
    _resolved: dict[str, Pmf | None] = field(
        default_factory=dict, init=False, compare=False, repr=False, hash=False
    )
    _generic_bases: set[str] = field(
        default_factory=set, init=False, compare=False, repr=False, hash=False
    )

    def __post_init__(self) -> None:
        if self.fallback is not None and not _has_error_mass(
            scale_error_odds(self.fallback.generic_pmf, 0.0, self.min_len)
        ):
            raise ValueError(
                "stutter_fallback generic pmf must keep error mass for runs of "
                f"{self.min_len} bases (an error delta >= -{self.min_len})"
            )

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Mapping[str, float]],
        min_len: int = 3,
        fallback: StutterFallback | None = None,
    ) -> StutterTable:
        pmfs: dict[str, Pmf] = {}
        for key, pmf in data.items():
            validate_stutter_key(key)
            pmfs[key] = parse_pmf(pmf, key, int(key.partition("|")[0][1:]))
        return cls(pmfs, min_len, fallback)

    def resolve(self, base: str, length: int, strand: str) -> Pmf | None:
        """The pmf used for a run; None when there is no entry and no fallback."""
        key = f"{base}{length}|{strand}"
        if key not in self._resolved:
            self._resolved[key] = self._lookup(base, length, strand)
        return self._resolved[key]

    def sample(self, base: str, length: int, strand: str, rng: random.Random) -> int:
        """Draw a length delta; 0 when the run resolves to no pmf."""
        pmf = self.resolve(base, length, strand)
        if not pmf:
            return 0
        draw = rng.random()
        cumulative = 0.0
        for delta, prob in pmf:
            cumulative += prob
            if draw < cumulative:
                return delta
        return pmf[-1][0]

    def _lookup(self, base: str, length: int, strand: str) -> Pmf | None:
        pmf = self.pmfs.get(f"{base}{length}|{strand}") or self.pmfs.get(f"{base}{length}|both")
        if pmf or self.fallback is None:
            return pmf
        key = f"{base}{length}|{strand}"
        fitted = self._fitted_lengths(base, strand)
        if fitted:
            shorter = [n for n in fitted if n <= length]
            longer = [n for n in fitted if n >= length]
            if shorter and longer:
                lo, hi = max(shorter), min(longer)
                between = interpolate_log_odds(
                    self.pmfs[fitted[lo]], lo, self.pmfs[fitted[hi]], hi, length
                )
                if between is not None:
                    logging.debug(
                        "Stutter %s: interpolated between %s and %s", key, fitted[lo], fitted[hi]
                    )
                    return between
            ref_len = max(shorter) if shorter else min(fitted)
            reference, source = self.pmfs[fitted[ref_len]], fitted[ref_len]
            logging.debug("Stutter %s has no fitted entry; extrapolating from %s", key, source)
        else:
            reference, ref_len = self.fallback.generic_pmf, self.fallback.generic_ref_len
            if base not in self._generic_bases:
                self._generic_bases.add(base)
                logging.warning(
                    "Stutter table has no fitted entry for base %s (first run: %s); runs of "
                    "%s use the profile's generic stutter pmf (ref_len %d)",
                    base,
                    key,
                    base,
                    ref_len,
                )
            logging.debug("Stutter %s uses the generic stutter pmf", key)
        shift = self.fallback.shift(length - ref_len)
        scaled = scale_error_odds(reference, shift, length)
        if _has_error_mass(scaled) or not fitted:
            return scaled
        # The reference's error deltas all remove more bases than this run has.
        logging.info("Stutter %s: reference has no usable error mass; using the generic pmf", key)
        generic_shift = self.fallback.shift(length - self.fallback.generic_ref_len)
        return scale_error_odds(self.fallback.generic_pmf, generic_shift, length)

    def _fitted_lengths(self, base: str, strand: str) -> dict[int, str]:
        """Run length -> key of fitted entries for ``base`` on ``strand`` (else ``both``)."""
        for wanted in (strand, "both"):
            found = {
                int(key.partition("|")[0][1:]): key
                for key in self.pmfs
                if key[0] == base and key.partition("|")[2] == wanted
            }
            if found:
                return found
        return {}

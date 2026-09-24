"""Strand-aware homopolymer stutter tables and their fallback for unfitted runs.

A :class:`StutterTable` maps keys such as ``"C7|+"`` (base and run length in
source orientation, read strand ``+``/``-`` or ``both``) to a probability mass
over length deltas. Fitted entries come from measured data.

Runs without a fitted entry are resolved by an optional
:class:`StutterFallback` (rule ``log_odds_linear``):

1. **Reference.** Use the fitted entries of the same base and strand (or of
   strand ``both`` when that base has no strand-specific entries). The
   reference is the longest fitted length that is not longer than the run; if
   every fitted length is longer, it is the shortest fitted length.
2. **Generic.** If the base has no fitted entries at all, use the profile's
   generic pmf at its ``ref_len`` and log a warning (once per base and table).
3. **Scaling.** The error odds ``P(delta != 0) / P(delta == 0)`` of the
   reference are multiplied by ``exp(log_odds_slope_per_base * (run length -
   reference length))``. The error shape (the pmf conditional on an error) is
   kept; deltas that would remove more bases than the run has are dropped.

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
FALLBACK_RULES = ("log_odds_linear",)
_FALLBACK_KEYS = {"rule", "log_odds_slope_per_base", "generic"}
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


@dataclass(frozen=True)
class StutterFallback:
    """Documented rule for runs without a fitted stutter entry (see module docstring)."""

    log_odds_slope_per_base: float
    generic_pmf: Pmf
    generic_ref_len: int
    rule: str = "log_odds_linear"

    def __post_init__(self) -> None:
        if self.rule not in FALLBACK_RULES:
            raise ValueError(f"stutter_fallback rule must be one of {FALLBACK_RULES}")
        slope = self.log_odds_slope_per_base
        if not math.isfinite(slope) or slope < 0:
            raise ValueError(f"log_odds_slope_per_base must be finite and >= 0, got {slope}")
        if self.generic_ref_len < 1:
            raise ValueError("stutter_fallback generic ref_len must be >= 1")
        if dict(self.generic_pmf).get(0, 0.0) >= 1.0:
            raise ValueError("stutter_fallback generic pmf must have error mass (P(0) < 1)")

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
        shift = self.fallback.log_odds_slope_per_base * (length - ref_len)
        return scale_error_odds(reference, shift, length)

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

"""Tests for fitting stutter tables from aggregate targets and validating them (#132)."""

from __future__ import annotations

import math
from collections import Counter

import pytest

from muc_one_up.read_simulator.stutter import StutterFallback, StutterTable
from muc_one_up.read_simulator.stutter_fit import (
    compare_to_targets,
    fit_log_odds_slope,
    fit_stutter_pmfs,
    pooled_generic_pmf,
    round_pmf,
    stutter_profile_section,
)


def _target(n: int, pmf: dict[int, float]) -> dict:
    return {
        "n": n,
        "p_correct": pmf.get(0, 0.0),
        "p_obs_minus_true": {str(d): p for d, p in pmf.items()},
    }


HP_TARGETS = {
    "C2|+": _target(9000, {-1: 0.02, 0: 0.98}),
    "C3|+": _target(5000, {-1: 0.05, 0: 0.9, 1: 0.05}),
    "C3|-": _target(5000, {-1: 0.04, 0: 0.96}),
    "C3|both": _target(10000, {-1: 0.045, 0: 0.93, 1: 0.025}),
    "C7|+": _target(3000, {-5: 0.01, -1: 0.39, 0: 0.5, 1: 0.1}),
    "C7|-": _target(3000, {-1: 0.1, 0: 0.9}),
    "C8|+": _target(854, {-7: 0.01, -1: 0.39, 0: 0.4, 1: 0.2}),
    "C8|-": _target(120, {-1: 0.2, 0: 0.8}),
    "G3|+": _target(4000, {-1: 0.02, 0: 0.95, 1: 0.03}),
}


class TestFitStutterPmfs:
    def test_selects_per_strand_keys_with_enough_observations(self) -> None:
        pmfs = fit_stutter_pmfs(HP_TARGETS, min_n=500, min_len=3, max_delta=6)
        assert set(pmfs) == {"C3|+", "C3|-", "C7|+", "C7|-", "C8|+", "G3|+"}

    def test_truncates_to_max_delta_and_renormalises(self) -> None:
        pmfs = fit_stutter_pmfs(HP_TARGETS, min_n=500, min_len=3, max_delta=6)
        assert -7 not in pmfs["C8|+"]
        assert pmfs["C8|+"][0] == pytest.approx(0.4 / 0.99)
        assert pmfs["C7|+"][-5] == pytest.approx(0.01)  # within +-6: kept
        assert sum(pmfs["C8|+"].values()) == pytest.approx(1.0)

    def test_fitted_pmfs_are_valid_stutter_entries(self) -> None:
        pmfs = fit_stutter_pmfs(HP_TARGETS, min_n=500, min_len=3, max_delta=3)
        StutterTable.from_dict({k: {str(d): p for d, p in v.items()} for k, v in pmfs.items()})


class TestLogOddsSlope:
    def test_median_least_squares_slope_over_groups_with_enough_lengths(self) -> None:
        pmfs = {
            "C3|+": {-1: 0.1, 0: 0.9},
            "C5|+": {-1: 0.5, 0: 0.5},
            "C7|+": {-1: 0.9, 0: 0.1},
            "G3|+": {-1: 0.1, 0: 0.9},  # single length: ignored
        }
        expected = (math.log(9.0) - math.log(1 / 9)) / 4
        assert fit_log_odds_slope(pmfs, min_lengths=3) == pytest.approx(expected)

    def test_negative_slopes_are_clamped_to_zero(self) -> None:
        pmfs = {"C3|+": {-1: 0.3, 0: 0.7}, "C4|+": {-1: 0.2, 0: 0.8}}
        assert fit_log_odds_slope(pmfs, min_lengths=2) == 0.0

    def test_requires_a_group_with_enough_lengths(self) -> None:
        with pytest.raises(ValueError, match="min_lengths"):
            fit_log_odds_slope({"C3|+": {-1: 0.3, 0: 0.7}}, min_lengths=2)


def test_generic_pmf_pools_fitted_keys_at_reference_length_by_n() -> None:
    pmfs = fit_stutter_pmfs(HP_TARGETS, min_n=500, min_len=3, max_delta=6)
    generic = pooled_generic_pmf(HP_TARGETS, pmfs, ref_len=3)
    # C3|+ (5000), C3|- (5000), G3|+ (4000)
    assert generic[0] == pytest.approx((5000 * 0.9 + 5000 * 0.96 + 4000 * 0.95) / 14000)
    assert sum(generic.values()) == pytest.approx(1.0)


def test_profile_section_is_loadable_and_records_provenance() -> None:
    molecules, provenance = stutter_profile_section(
        HP_TARGETS,
        min_n=500,
        min_len=3,
        max_delta=6,
        min_lengths_for_slope=2,
        generic_ref_len=3,
        max_extrapolation_bases=3,
    )
    fallback = StutterFallback.from_dict(molecules["stutter_fallback"])
    table = StutterTable.from_dict(molecules["stutter"], fallback=fallback)
    assert table.resolve("A", 4, "+") is not None
    assert provenance["min_target_n"] == 500
    assert molecules["stutter_fallback"]["max_extrapolation_bases"] == 3
    assert provenance["fitted_keys"] == sorted(molecules["stutter"])
    assert "C8|-" in provenance["excluded_keys"]
    for pmf in molecules["stutter"].values():
        assert sum(pmf.values()) == pytest.approx(1.0, abs=1e-9)


class TestCompareToTargets:
    def test_reports_every_per_strand_target_key(self) -> None:
        observed = {
            "C7|+": Counter({0: 50, -1: 40, 1: 10}),
            "C8|+": Counter({0: 20, -1: 60, 1: 20}),
        }
        report = compare_to_targets(observed, HP_TARGETS, max_delta=6, tolerance=0.03)
        assert set(report) == {k for k in HP_TARGETS if not k.endswith("|both")}
        c7 = report["C7|+"]
        assert c7["n_sim"] == 100 and c7["p_correct_sim"] == pytest.approx(0.5)
        assert c7["p_correct_target"] == pytest.approx(0.5)
        assert c7["max_abs_residual"] == pytest.approx(0.01)
        assert c7["within_tolerance"] is True
        assert report["C8|+"]["within_tolerance"] is False  # p_correct 0.2 vs 0.4
        assert report["G3|+"]["status"] == "not_in_template"
        assert report["G3|+"]["within_tolerance"] is None


class TestNumericalGuards:
    def test_target_without_mass_in_range_is_an_error_naming_the_key(self) -> None:
        targets = {"C8|+": _target(900, {-7: 1.0})}
        with pytest.raises(ValueError, match="C8"):
            fit_stutter_pmfs(targets, min_n=500, min_len=3, max_delta=6)

    def test_rounding_never_makes_p0_negative(self) -> None:
        pmf = {0: 0.000002, -1: 0.333336, 1: 0.333336, 2: 0.333326}  # rounds to 1.00001
        rounded = round_pmf(pmf)
        assert all(p >= 0 for p in rounded.values())
        assert sum(rounded.values()) == pytest.approx(1.0, abs=1e-9)

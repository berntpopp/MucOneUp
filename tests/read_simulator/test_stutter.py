"""Tests for homopolymer stutter tables and their fallback for unfitted runs (#132)."""

from __future__ import annotations

import logging
import math
import random

import pytest

from muc_one_up.read_simulator.molecules import MoleculeModel
from muc_one_up.read_simulator.stutter import StutterFallback, StutterTable, scale_error_odds

SLOPE = math.log(2.0)  # error odds double per extra base
FALLBACK = {
    "rule": "log_odds_linear",
    "log_odds_slope_per_base": SLOPE,
    "generic": {"ref_len": 3, "pmf": {"-1": 0.04, "0": 0.9, "1": 0.06}},
}
TABLE = {
    "C3|+": {"-1": 0.1, "0": 0.9},
    "C7|+": {"-2": 0.1, "-1": 0.3, "0": 0.5, "1": 0.1},
    "C7|-": {"-1": 0.2, "0": 0.8},
    "G4|both": {"-1": 0.25, "0": 0.75},
}


def _table() -> StutterTable:
    return StutterTable.from_dict(TABLE, fallback=StutterFallback.from_dict(FALLBACK))


def _p0(pmf: tuple[tuple[int, float], ...]) -> float:
    return dict(pmf).get(0, 0.0)


class TestWithoutFallback:
    def test_missing_key_draws_no_delta(self) -> None:
        """Tables without a fallback keep the pbsim3-path behaviour: no extra stutter."""
        table = StutterTable.from_dict({"C7|+": {"-1": 1.0}})
        assert table.sample("G", 4, "+", random.Random(1)) == 0
        assert table.resolve("G", 4, "+") is None


class TestFallbackResolution:
    def test_fitted_key_is_used_unchanged(self) -> None:
        assert _table().resolve("C", 7, "+") == ((-2, 0.1), (-1, 0.3), (0, 0.5), (1, 0.1))

    def test_longer_run_extrapolates_from_longest_fitted_length(self) -> None:
        pmf = _table().resolve("C", 8, "+")
        assert pmf is not None
        odds = 0.5 / 0.5 * 2.0  # C7|+ error odds, doubled for one extra base
        assert _p0(pmf) == pytest.approx(1 / (1 + odds))
        errors = {d: p for d, p in pmf if d != 0}
        # the error shape (conditional on an error) is kept from C7|+
        assert errors[-1] / errors[-2] == pytest.approx(3.0)
        assert errors[-1] / errors[1] == pytest.approx(3.0)
        assert sum(p for _, p in pmf) == pytest.approx(1.0)

    def test_strand_is_kept_when_extrapolating(self) -> None:
        pmf = _table().resolve("C", 9, "-")
        assert pmf is not None
        assert _p0(pmf) == pytest.approx(1 / (1 + 0.25 * 4.0))

    def test_gap_uses_longest_fitted_length_below(self) -> None:
        pmf = _table().resolve("C", 5, "+")  # fitted: C3 and C7
        assert pmf is not None
        assert _p0(pmf) == pytest.approx(1 / (1 + (0.1 / 0.9) * 4.0))

    def test_shorter_than_fitted_scales_down_from_shortest(self) -> None:
        pmf = _table().resolve("C", 3, "-")  # only C7|- fitted
        assert pmf is not None
        assert _p0(pmf) == pytest.approx(1 / (1 + 0.25 / 16.0))

    def test_both_strand_entries_serve_as_reference(self) -> None:
        pmf = _table().resolve("G", 5, "+")
        assert pmf is not None
        assert _p0(pmf) == pytest.approx(1 / (1 + (0.25 / 0.75) * 2.0))

    def test_unknown_base_uses_generic_pmf_and_warns_once_per_base(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        table = _table()
        with caplog.at_level(logging.WARNING):
            pmf = table.resolve("A", 4, "+")
            table.resolve("A", 4, "+")
            table.resolve("A", 6, "-")
            table.resolve("T", 5, "+")
        assert pmf is not None
        assert _p0(pmf) == pytest.approx(1 / (1 + (0.1 / 0.9) * 2.0))
        warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 2
        assert "base A" in warnings[0] and "A4|+" in warnings[0] and "generic" in warnings[0]
        assert "base T" in warnings[1]

    def test_deltas_removing_more_than_the_run_are_dropped(self) -> None:
        table = StutterTable.from_dict(
            {"T6|+": {"-4": 0.1, "-1": 0.1, "0": 0.8}},
            fallback=StutterFallback.from_dict({**FALLBACK, "log_odds_slope_per_base": 0.0}),
        )
        pmf = table.resolve("T", 3, "+")
        assert pmf is not None
        assert min(d for d, _ in pmf) >= -3
        assert _p0(pmf) == pytest.approx(0.8 / 0.9)

    def test_every_protected_run_gets_an_error_mass(self) -> None:
        table = _table()
        for base in "ACGT":
            for length in range(3, 15):
                for strand in "+-":
                    pmf = table.resolve(base, length, strand)
                    assert pmf is not None and _p0(pmf) < 1.0, f"{base}{length}|{strand}"

    def test_sampling_follows_the_extrapolated_pmf(self) -> None:
        table = _table()
        rng = random.Random(11)
        draws = [table.sample("C", 8, "+", rng) for _ in range(20000)]
        assert draws.count(0) / len(draws) == pytest.approx(1 / 3, abs=0.015)


class TestFallbackValidation:
    @pytest.mark.parametrize(
        "mutation,match",
        [
            ({"rule": "nearest"}, "rule"),
            ({"log_odds_slope_per_base": -0.1}, "slope"),
            ({"log_odds_slope_per_base": float("nan")}, "slope"),
            ({"generic": {"ref_len": 3, "pmf": {"0": 1.0}}}, "error"),
            ({"generic": {"ref_len": 3, "pmf": {"0": 0.5, "1": 0.2}}}, "sum"),
            ({"generic": {"ref_len": 2, "pmf": {"-3": 0.1, "0": 0.9}}}, "remove"),
            ({"generic": {"ref_len": 0, "pmf": {"-1": 0.1, "0": 0.9}}}, "ref_len"),
            ({"surprise": 1}, "unknown"),
        ],
    )
    def test_invalid_fallback_rejected(self, mutation: dict, match: str) -> None:
        with pytest.raises(ValueError, match=match):
            StutterFallback.from_dict({**FALLBACK, **mutation})

    def test_molecule_model_reads_fallback_next_to_stutter(self) -> None:
        model = MoleculeModel.from_dict({"stutter": TABLE, "stutter_fallback": FALLBACK})
        assert model.stutter is not None and model.stutter.fallback is not None
        assert model.stutter.resolve("C", 8, "+") is not None

    def test_fallback_without_table_rejected(self) -> None:
        with pytest.raises(ValueError, match="stutter_fallback"):
            MoleculeModel.from_dict({"stutter_fallback": FALLBACK})


class TestNoErrorFreeFallback:
    """A resolved fallback pmf always keeps error mass (#132 review)."""

    def test_reference_losing_all_error_deltas_uses_generic(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        table = StutterTable.from_dict(
            {"C7|+": {"-6": 0.3, "0": 0.7}}, fallback=StutterFallback.from_dict(FALLBACK)
        )
        for length in (3, 4, 5):
            with caplog.at_level(logging.INFO):
                pmf = table.resolve("C", length, "+")
            assert pmf is not None and _p0(pmf) < 1.0, length
            assert min(d for d, _ in pmf) >= -length
        assert "generic" in caplog.text

    def test_generic_pmf_must_keep_error_mass_at_min_len(self) -> None:
        fallback = StutterFallback.from_dict(
            {**FALLBACK, "generic": {"ref_len": 5, "pmf": {"-5": 0.1, "0": 0.9}}}
        )
        with pytest.raises(ValueError, match="generic"):
            StutterTable.from_dict({"C7|+": {"0": 0.5, "-1": 0.5}}, fallback=fallback)


class TestScaleErrorOdds:
    def test_error_free_reference_has_no_error_mass_to_scale(self) -> None:
        """The pure scaler cannot invent errors; StutterTable then uses the generic pmf."""
        assert scale_error_odds(((0, 1.0),), 2.0, 5) == ((0, 1.0),)

    def test_error_only_reference_stays_error_only(self) -> None:
        assert scale_error_odds(((-1, 0.5), (1, 0.5)), 2.0, 5) == ((-1, 0.5), (1, 0.5))

    def test_zero_shift_is_identity_up_to_dropped_deltas(self) -> None:
        pmf = ((-1, 0.2), (0, 0.7), (1, 0.1))
        scaled = scale_error_odds(pmf, 0.0, 5)
        assert [d for d, _ in scaled] == [-1, 0, 1]
        assert [p for _, p in scaled] == pytest.approx([0.2, 0.7, 0.1])

    def test_malformed_generic_block_rejected(self) -> None:
        with pytest.raises(ValueError, match="generic"):
            StutterFallback.from_dict({**FALLBACK, "generic": {"pmf": {"0": 1.0}}})

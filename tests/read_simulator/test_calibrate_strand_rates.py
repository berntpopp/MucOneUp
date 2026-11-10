"""Convergence check for helpers/calibrate_read_profile.py --fit-strand-rates (#131)."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest

pytest.importorskip("edlib")

from muc_one_up.read_simulator.molecules import MoleculeModel

HELPER = Path(__file__).resolve().parents[2] / "helpers" / "calibrate_read_profile.py"


def _helper():
    spec = importlib.util.spec_from_file_location("calibrate_read_profile", HELPER)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fit_strand_rates_moves_towards_asymmetric_targets() -> None:
    errors = {
        "mismatch_rate": 0.01,
        "insertion_rate": 0.01,
        "deletion_rate": 0.01,
        "insertion_len_pmf": {"1": 1.0},
        "deletion_len_pmf": {"1": 1.0},
        "read_error_sigma": 0.0,
    }
    target = {"+": (0.02, 0.005, 0.015), "-": (0.004, 0.012, 0.002)}
    targets = {
        "error_rates_per_ref_base": {
            f"{s}:{kind}": {"median": value}
            for s, values in target.items()
            for kind, value in zip(("mismatch", "ins", "del"), values, strict=True)
        }
    }
    rng = random.Random(0)
    template = "".join(rng.choice("ACGT") for _ in range(2000))
    fit = _helper().fit_strand_rates(
        errors, targets, template, MoleculeModel(forward_frac=0.5), 200, 1, 8, 2
    )

    def gap(sim: dict[str, float]) -> float:
        return sum(
            abs(sim[k] - t["median"]) for k, t in targets["error_rates_per_ref_base"].items()
        )

    rounds = fit["report"]["simulated_per_round"]
    assert len(rounds) == 3
    assert gap(rounds[-1]) < 0.25 * gap(rounds[0])
    rates = fit["strand_rates"]
    assert rates["+"]["mismatch_rate"] > rates["-"]["mismatch_rate"]
    assert rates["+"]["deletion_rate"] > rates["-"]["deletion_rate"]


def test_refit_unit_stutter_converges_on_the_profiler_statistic() -> None:
    rng = random.Random(0)
    units = 20
    vntr = "".join(rng.choice("ACGT") for _ in range(60 * units))  # distinct units: unique anchors
    template = "ACGTACGTAC" + vntr + "GTCAGTCAGT"
    target = {"minus2_per_onpeak": 0.02, "minus1_per_onpeak": 0.08, "plus1_per_onpeak": 0.01}
    targets = {
        "per_allele_stutter_and_length_noise": [
            {"lib": "x", "n_on_peak": 1000, "units": units, **target}
        ]
    }
    errors = {
        "mismatch_rate": 0.005,
        "insertion_rate": 0.005,
        "deletion_rate": 0.005,
        "insertion_len_pmf": {"1": 1.0},
        "deletion_len_pmf": {"1": 1.0},
        "read_error_sigma": 0.0,
    }
    helper = _helper()
    run = helper.EmpiricalSequencer(helper.EmpiricalErrorModel.from_dict(errors))
    start = {"unit_stutter_rate": 0.01, "unit_stutter_pmf": {"-2": 0.2, "-1": 0.7, "1": 0.1}}
    section, record = helper.refit_unit_stutter(
        start, {"forward_frac": 0.5}, targets, template, (10, 10 + 60 * units), run, 1500, 1, 2
    )
    rounds = record["rounds"]
    assert len(rounds) == 3 and rounds[-1]["unit_stutter_rate"] == section["unit_stutter_rate"]

    def gap(r: dict) -> float:
        return abs(r["per_onpeak"]["-1"] - 0.08) + abs(r["per_onpeak"]["-2"] - 0.02)

    assert gap(rounds[-1]) < 0.3 * gap(rounds[0])
    assert section["unit_stutter_rate"] < start["unit_stutter_rate"]

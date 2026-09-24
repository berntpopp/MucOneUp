"""No homopolymer run in a simulated haplotype is error-free by construction (#132).

The empirical error channel protects every run >= ``hp_min_len`` from base-level
errors, so each such run must receive stutter. Before #132 the dupC C8 run and
all A/T runs of the ONT profiles had no stutter entry and were read perfectly.
"""

from __future__ import annotations

import random

import pytest

from muc_one_up.read_simulator.molecules import (
    MoleculeModel,
    build_amplicon_molecules,
    homopolymer_runs,
)
from muc_one_up.read_simulator.read_profiles import list_builtin_profiles, load_read_profile

UNIT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"  # X: C7 at 52-58
DUPC_UNIT = UNIT[:52] + "C" + UNIT[52:]  # dupC: C8 at 52-59
A4_LINKER = "GTCTAAAAGTCTGA"
HAPLOTYPE = UNIT + DUPC_UNIT + A4_LINKER + UNIT
C8_START = len(UNIT) + 52
A4_START = len(UNIT) + len(DUPC_UNIT) + A4_LINKER.index("AAAA")

# PRJEB92208 amplicon aggregates (hp_P_obs_given_true, p_correct): C8|+ n=854, C8|- n=1144
TARGET_C8_P_CORRECT = {"+": 0.362, "-": 0.747}
N_MOLECULES = 8000
SEED = 132
TOLERANCE = 0.03  # >= 4 binomial standard errors at ~4000 molecules per strand

EMPIRICAL_PROFILES = [
    name for name in list_builtin_profiles() if load_read_profile(name).errors is not None
]


def _unedited_fraction(start: int, stutter_model: MoleculeModel) -> dict[str, float]:
    molecules = build_amplicon_molecules(
        [HAPLOTYPE], [N_MOLECULES], stutter_model, random.Random(SEED)
    )
    total = dict.fromkeys("+-", 0)
    unedited = dict.fromkeys("+-", 0)
    for molecule in molecules:
        total[molecule.strand] += 1
        if not any(edit.pos == start for edit in molecule.hp_edits):
            unedited[molecule.strand] += 1
    return {strand: unedited[strand] / total[strand] for strand in "+-"}


def test_haplotype_contains_the_runs_under_test() -> None:
    runs = {(s, e - s, b) for s, e, b in homopolymer_runs(HAPLOTYPE, 3)}
    assert (C8_START, 8, "C") in runs
    assert (A4_START, 4, "A") in runs


@pytest.mark.parametrize("name", EMPIRICAL_PROFILES)
def test_no_protected_run_is_error_free_by_construction(name: str) -> None:
    profile = load_read_profile(name)
    table = profile.molecules.stutter
    assert profile.errors is not None and table is not None
    for _start, end, base in homopolymer_runs(HAPLOTYPE, profile.errors.hp_min_len):
        for strand in "+-":
            pmf = table.resolve(base, end - _start, strand)
            assert pmf is not None, f"{base}{end - _start}|{strand} has no stutter pmf"
            assert dict(pmf).get(0, 0.0) < 1.0, f"{base}{end - _start}|{strand} is error-free"


def test_dupc_c8_is_stuttered_at_the_target_rate() -> None:
    profile = load_read_profile("ont_r10_sup_amplicon_v1")
    model = MoleculeModel(forward_frac=0.5, stutter=profile.molecules.stutter)
    p_correct = _unedited_fraction(C8_START, model)
    for strand, target in TARGET_C8_P_CORRECT.items():
        assert p_correct[strand] == pytest.approx(target, abs=TOLERANCE), strand


@pytest.mark.parametrize("name", EMPIRICAL_PROFILES)
def test_a4_run_receives_stutter(name: str) -> None:
    profile = load_read_profile(name)
    model = MoleculeModel(forward_frac=0.5, stutter=profile.molecules.stutter)
    p_correct = _unedited_fraction(A4_START, model)
    for strand in "+-":
        assert p_correct[strand] < 1.0, strand

"""Tests for the template molecule model (strand, PCR artefacts, homopolymer stutter)."""

from __future__ import annotations

import random

import pytest

from muc_one_up.read_simulator.molecules import (
    MoleculeModel,
    StutterTable,
    apply_stutter,
    build_amplicon_molecules,
    build_fragment_molecules,
    reverse_complement,
)

UNIT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"  # X: C7 at 52-58
AMP_A = UNIT * 30
AMP_B = UNIT * 45


def test_reverse_complement() -> None:
    assert reverse_complement("AACGTN") == "NACGTT"


class TestStutterTable:
    def test_rejects_pmf_not_summing_to_one(self) -> None:
        with pytest.raises(ValueError, match="sum"):
            StutterTable.from_dict({"C7|+": {"0": 0.5, "1": 0.2}})

    def test_rejects_malformed_key(self) -> None:
        with pytest.raises(ValueError, match="key"):
            StutterTable.from_dict({"C7": {"0": 1.0}})

    def test_missing_key_means_no_stutter(self) -> None:
        table = StutterTable.from_dict({"C7|+": {"-1": 1.0}})
        assert table.sample("G", 4, "+", random.Random(1)) == 0


class TestApplyStutter:
    def test_certain_deletion_edits_every_matching_run(self) -> None:
        table = StutterTable.from_dict({"C7|+": {"-1": 1.0}})
        seq, edits = apply_stutter(UNIT * 3, table, "+", random.Random(1))
        assert len(seq) == len(UNIT) * 3 - 3
        assert [(e.base, e.true_len, e.new_len) for e in edits] == [("C", 7, 6)] * 3
        assert [e.pos for e in edits] == [52, 112, 172]

    def test_strand_selects_table(self) -> None:
        table = StutterTable.from_dict({"C7|-": {"1": 1.0}})
        assert apply_stutter(UNIT, table, "+", random.Random(1))[0] == UNIT
        assert "C" * 8 in apply_stutter(UNIT, table, "-", random.Random(1))[0]

    def test_empirical_rate_matches_pmf(self) -> None:
        table = StutterTable.from_dict({"C7|+": {"-1": 0.3, "0": 0.6, "1": 0.1}})
        rng = random.Random(7)
        deltas = [apply_stutter(UNIT, table, "+", rng)[1] for _ in range(4000)]
        minus = sum(1 for d in deltas if d and d[0].new_len == 6) / 4000
        plus = sum(1 for d in deltas if d and d[0].new_len == 8) / 4000
        assert minus == pytest.approx(0.3, abs=0.03)
        assert plus == pytest.approx(0.1, abs=0.02)


class TestMoleculeModel:
    def test_defaults_are_a_no_op(self) -> None:
        mols = build_amplicon_molecules([AMP_A], [20], MoleculeModel(), random.Random(1))
        assert all(m.seq == AMP_A and m.kind == "full" and m.strand == "+" for m in mols)

    def test_from_dict_rejects_unknown_keys(self) -> None:
        with pytest.raises(ValueError, match="unknown"):
            MoleculeModel.from_dict({"smear": 0.1})

    @pytest.mark.parametrize("key,value", [("forward_frac", 1.5), ("smear_rate", -0.1)])
    def test_from_dict_validates_ranges(self, key: str, value: float) -> None:
        with pytest.raises(ValueError, match=key):
            MoleculeModel.from_dict({key: value})


class TestAmpliconMolecules:
    def test_counts_ids_and_haplotypes(self) -> None:
        mols = build_amplicon_molecules([AMP_A, AMP_B], [30, 10], MoleculeModel(), random.Random(1))
        assert [m.id for m in mols] == list(range(1, 41))
        assert sum(m.hap == 1 for m in mols) == 30 and sum(m.hap == 2 for m in mols) == 10

    def test_strand_mix(self) -> None:
        model = MoleculeModel(forward_frac=0.5)
        mols = build_amplicon_molecules([AMP_A], [2000], model, random.Random(3))
        forward = sum(m.strand == "+" for m in mols) / len(mols)
        assert forward == pytest.approx(0.5, abs=0.05)
        assert all(m.seq == reverse_complement(AMP_A) for m in mols if m.strand == "-")

    def test_smear_products_are_shorter_single_junction_deletions(self) -> None:
        model = MoleculeModel(smear_rate=0.3, smear_min_keep=0.2)
        mols = build_amplicon_molecules([AMP_B], [1000], model, random.Random(5))
        smear = [m for m in mols if m.kind == "smear"]
        assert len(smear) / len(mols) == pytest.approx(0.3, abs=0.05)
        for m in smear:
            start, resume = (int(x) for x in m.detail.removeprefix("deletion:").split("-"))
            assert m.seq == AMP_B[:start] + AMP_B[resume:]
            assert 0.2 * len(AMP_B) - 1 <= len(m.seq) < len(AMP_B)

    def test_chimeras_join_two_haplotypes_at_homologous_fraction(self) -> None:
        model = MoleculeModel(chimera_rate=1.0)
        mols = build_amplicon_molecules([AMP_A, AMP_B], [50, 50], model, random.Random(2))
        chim = [m for m in mols if m.kind == "chimera"]
        assert chim and all(len(AMP_A) <= len(m.seq) <= len(AMP_B) for m in chim)

    def test_no_chimeras_for_single_haplotype(self) -> None:
        model = MoleculeModel(chimera_rate=1.0)
        mols = build_amplicon_molecules([AMP_A], [50], model, random.Random(2))
        assert all(m.kind == "full" for m in mols)

    def test_concatemers_and_offtarget(self) -> None:
        model = MoleculeModel(concatemer_rate=0.1, offtarget_frac=0.25, offtarget_median_bp=300)
        mols = build_amplicon_molecules(
            [AMP_A], [400], model, random.Random(4), offtarget_source="ACGT" * 2000
        )
        assert sum(m.kind == "concatemer" for m in mols) == pytest.approx(40, abs=15)
        off = [m for m in mols if m.kind == "offtarget"]
        assert len(off) == round(400 * 0.25 / 0.75)  # offtarget_frac of the final total
        assert all(m.hap == 0 for m in off)

    def test_offtarget_requires_source(self) -> None:
        with pytest.raises(ValueError, match="offtarget_source"):
            build_amplicon_molecules(
                [AMP_A], [10], MoleculeModel(offtarget_frac=0.2), random.Random(1)
            )

    def test_deterministic_for_seed(self) -> None:
        model = MoleculeModel(forward_frac=0.5, smear_rate=0.2, chimera_rate=0.05)
        a = build_amplicon_molecules([AMP_A, AMP_B], [100, 60], model, random.Random(9))
        b = build_amplicon_molecules([AMP_A, AMP_B], [100, 60], model, random.Random(9))
        assert a == b


class TestFragmentMolecules:
    def test_lengths_haplotype_balance_and_bounds(self) -> None:
        src = ["A" * 20000, "C" * 20000]
        mols = build_fragment_molecules(src, 2000, 5000, 0.4, MoleculeModel(), random.Random(1))
        assert len(mols) == 2000
        assert sum(m.hap == 1 for m in mols) / 2000 == pytest.approx(0.5, abs=0.05)
        lengths = sorted(len(m.seq) for m in mols)
        assert 3500 < lengths[1000] <= 5200
        assert all(0 <= m.src_start < m.src_end <= 20000 and m.kind == "fragment" for m in mols)

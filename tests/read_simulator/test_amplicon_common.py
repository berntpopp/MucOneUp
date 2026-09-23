"""Tests for shared amplicon extraction and preparation."""

from pathlib import Path

import pytest

from muc_one_up.read_simulator.amplicon_common import AmpliconPrep


@pytest.fixture
def muc1_primers():
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def diploid_fasta_with_primers(tmp_path, muc1_primers):
    """Create a diploid FASTA with primer sites in both haplotypes."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())

    vntr1 = "ACGT" * 400
    seq1 = "N" * 50 + fwd + vntr1 + rev_rc + "N" * 50

    vntr2 = "ACGT" * 800
    seq2 = "N" * 50 + fwd + vntr2 + rev_rc + "N" * 50

    fasta = tmp_path / "diploid.fa"
    fasta.write_text(f">hap1\n{seq1}\n>hap2\n{seq2}\n")
    return fasta


class TestAmpliconPrep:
    def test_dataclass_fields(self):
        """AmpliconPrep has required fields."""
        prep = AmpliconPrep(
            allele_templates=[Path("/a.fa")],
            allele_coverages=[100],
            output_dir=Path("/out"),
            output_base="test",
            intermediate_files=[],
            is_diploid=False,
        )
        assert prep.is_diploid is False
        assert len(prep.allele_templates) == 1


class TestExtractAndPrepare:
    def test_diploid_returns_two_templates(
        self, diploid_fasta_with_primers, tmp_path, muc1_primers
    ):
        """Diploid input produces 2 template FASTAs with PCR bias split."""
        from muc_one_up.read_simulator.amplicon_common import extract_and_prepare_amplicons

        prep = extract_and_prepare_amplicons(
            input_fa=str(diploid_fasta_with_primers),
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            total_coverage=100,
            work_dir=tmp_path,
            expected_product_range=None,
            pcr_bias_config={},
            seed=42,
        )

        assert prep.is_diploid is True
        assert len(prep.allele_templates) == 2
        assert len(prep.allele_coverages) == 2
        assert sum(prep.allele_coverages) >= 2  # each allele gets at least 1
        for t in prep.allele_templates:
            assert t.exists()


class TestTruthTrackedHelpers:
    """Shared truth-tracked helpers used by both amplicon pipelines (#100, #103)."""

    def test_model_selection(self, tmp_path):
        import json

        from muc_one_up.read_simulator.amplicon_common import truth_tracked_model

        assert truth_tracked_model({}, tracking_requested=False) is None
        assert truth_tracked_model({}, tracking_requested=True).forward_frac == 1.0
        path = tmp_path / "p.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "name": "p",
                    "platform": "ont",
                    "molecules": {"smear_rate": 0.2},
                }
            )
        )
        model = truth_tracked_model(
            {"read_model": {"profile": str(path)}}, tracking_requested=False
        )
        assert model.smear_rate == 0.2

    def test_offtarget_source_excludes_amplicon(self, tmp_path):
        from unittest.mock import patch

        from muc_one_up.read_simulator.amplicon_common import (
            AmpliconPrep,
            simulate_truth_tracked_amplicons,
        )
        from muc_one_up.read_simulator.molecule_pipeline import PbsimRun
        from muc_one_up.read_simulator.molecules import MoleculeModel

        prep = AmpliconPrep(
            allele_templates=[],
            allele_coverages=[5],
            output_dir=tmp_path,
            output_base="x",
            amplicon_sequences=["CCCCGGGG"],
            haplotype_sequences=["AAAACCCCGGGGTTTT"],
        )
        run = PbsimRun("pbsim", "samtools", "qshmm", "m.model")
        model = MoleculeModel(offtarget_frac=0.5, offtarget_median_bp=3)
        with patch("muc_one_up.read_simulator.amplicon_common.simulate_molecule_reads") as sim:
            sim.return_value = 10
            simulate_truth_tracked_amplicons(
                prep, model, run, tmp_path, tmp_path / "o.fq", tmp_path / "t.tsv.gz", "x", 1
            )
        molecules = sim.call_args.args[0]
        off = [m for m in molecules if m.kind == "offtarget"]
        assert len(off) == 5
        # No piece joins the left and right flank across the deleted amplicon (#112).
        assert all(set(m.seq) in ({"A"}, {"T"}) for m in off)
        for m in off:
            assert m.hap == 1
            span = (0, 4) if set(m.seq) == {"A"} else (12, 16)
            assert span[0] <= m.src_start < m.src_end <= span[1]

    def test_offtarget_without_flanks_raises_clear_error(self, tmp_path):
        from muc_one_up.read_simulator.amplicon_common import offtarget_intervals

        with pytest.raises(ValueError, match="no flanking sequence"):
            offtarget_intervals(["CCCCGGGG"], ["CCCCGGGG"])

    def test_offtarget_intervals_per_haplotype(self):
        from muc_one_up.read_simulator.amplicon_common import offtarget_intervals

        got = offtarget_intervals(["AACCGGT", "TTCCGG"], ["CCGG", "CCGG"])
        assert [(i.hap, i.start, i.seq) for i in got] == [(1, 0, "AA"), (1, 6, "T"), (2, 0, "TT")]

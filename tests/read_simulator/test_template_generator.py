"""Tests for template FASTA generator."""

import pytest
from Bio import SeqIO

from muc_one_up.read_simulator.utils.template_generator import generate_template_fasta


@pytest.fixture
def amplicon_fasta(tmp_path):
    """Create a single-sequence amplicon FASTA."""
    fasta = tmp_path / "amplicon.fa"
    fasta.write_text(">hap1_amplicon\nACGTACGTACGTACGT\n")
    return fasta


class TestGenerateTemplateFasta:
    """Tests for generate_template_fasta()."""

    def test_correct_record_count(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 100, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        assert len(records) == 100

    def test_all_records_have_same_sequence(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 10, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        seqs = [str(r.seq) for r in records]
        assert all(s == seqs[0] for s in seqs)

    def test_record_naming(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 3, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        assert records[0].id == "amplicon_copy_001"
        assert records[1].id == "amplicon_copy_002"
        assert records[2].id == "amplicon_copy_003"

    def test_single_copy(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 1, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        assert len(records) == 1

    def test_returns_output_path(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        result = generate_template_fasta(str(amplicon_fasta), 5, str(output))
        assert result == str(output)

    def test_creates_parent_directory(self, amplicon_fasta, tmp_path):
        output = tmp_path / "subdir" / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 5, str(output))
        assert output.exists()

    def test_zero_copies_raises(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        with pytest.raises(ValueError, match="num_copies must be"):
            generate_template_fasta(str(amplicon_fasta), 0, str(output))

    def test_negative_copies_raises(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        with pytest.raises(ValueError, match="num_copies must be"):
            generate_template_fasta(str(amplicon_fasta), -1, str(output))

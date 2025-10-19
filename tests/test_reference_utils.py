"""Comprehensive tests for reference_utils module.

Tests cover:
- Reference information extraction
- Diploid detection
- Haplotype extraction
- Edge cases and error handling
"""

import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from muc_one_up.exceptions import ValidationError
from muc_one_up.read_simulator.utils.reference_utils import (
    extract_haplotypes,
    get_reference_info,
    is_diploid_reference,
    validate_reference_compatibility,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def haploid_fasta(temp_dir):
    """Create a haploid FASTA file (1 sequence)."""
    fasta_path = temp_dir / "haploid.fa"
    record = SeqRecord(Seq("ATCGATCGATCG"), id="haplotype1", description="")
    SeqIO.write([record], fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def diploid_fasta(temp_dir):
    """Create a diploid FASTA file (2 sequences)."""
    fasta_path = temp_dir / "diploid.fa"
    records = [
        SeqRecord(Seq("ATCGATCGATCG"), id="haplotype_1", description=""),
        SeqRecord(Seq("GCTAGCTAGCTA"), id="haplotype_2", description=""),
    ]
    SeqIO.write(records, fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def diploid_unequal_fasta(temp_dir):
    """Create a diploid FASTA with unequal haplotype lengths."""
    fasta_path = temp_dir / "diploid_unequal.fa"
    records = [
        SeqRecord(Seq("ATCG" * 5), id="short_haplotype", description=""),  # 20 bp
        SeqRecord(Seq("GCTA" * 25), id="long_haplotype", description=""),  # 100 bp
    ]
    SeqIO.write(records, fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def triploid_fasta(temp_dir):
    """Create a FASTA file with 3 sequences."""
    fasta_path = temp_dir / "triploid.fa"
    records = [
        SeqRecord(Seq("AAAA"), id="seq1", description=""),
        SeqRecord(Seq("TTTT"), id="seq2", description=""),
        SeqRecord(Seq("GGGG"), id="seq3", description=""),
    ]
    SeqIO.write(records, fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def empty_fasta(temp_dir):
    """Create an empty FASTA file."""
    fasta_path = temp_dir / "empty.fa"
    fasta_path.write_text("")
    return fasta_path


# Tests for get_reference_info()


def test_get_reference_info_haploid(haploid_fasta):
    """Test reference info extraction for haploid reference."""
    info = get_reference_info(haploid_fasta)

    assert info.num_sequences == 1
    assert info.sequence_ids == ["haplotype1"]
    assert info.sequence_lengths == [12]
    assert info.total_length == 12
    assert info.is_diploid is False


def test_get_reference_info_diploid(diploid_fasta):
    """Test reference info extraction for diploid reference."""
    info = get_reference_info(diploid_fasta)

    assert info.num_sequences == 2
    assert info.sequence_ids == ["haplotype_1", "haplotype_2"]
    assert info.sequence_lengths == [12, 12]
    assert info.total_length == 24
    assert info.is_diploid is True


def test_get_reference_info_diploid_unequal(diploid_unequal_fasta):
    """Test reference info for diploid with unequal lengths."""
    info = get_reference_info(diploid_unequal_fasta)

    assert info.num_sequences == 2
    assert info.sequence_ids == ["short_haplotype", "long_haplotype"]
    assert info.sequence_lengths == [20, 100]
    assert info.total_length == 120
    assert info.is_diploid is True


def test_get_reference_info_triploid(triploid_fasta):
    """Test reference info for file with 3 sequences."""
    info = get_reference_info(triploid_fasta)

    assert info.num_sequences == 3
    assert info.sequence_ids == ["seq1", "seq2", "seq3"]
    assert info.sequence_lengths == [4, 4, 4]
    assert info.total_length == 12
    assert info.is_diploid is False


def test_get_reference_info_missing_file(temp_dir):
    """Test that missing file raises ValidationError."""
    missing_file = temp_dir / "nonexistent.fa"

    with pytest.raises(ValidationError, match="not found"):
        get_reference_info(missing_file)


def test_get_reference_info_empty_file(empty_fasta):
    """Test that empty FASTA raises ValidationError."""
    with pytest.raises(ValidationError, match="empty"):
        get_reference_info(empty_fasta)


def test_get_reference_info_directory(temp_dir):
    """Test that directory path raises ValidationError."""
    with pytest.raises(ValidationError, match="not a file"):
        get_reference_info(temp_dir)


def test_get_reference_info_accepts_string_path(diploid_fasta):
    """Test that function accepts string path (not just Path)."""
    info = get_reference_info(str(diploid_fasta))
    assert info.is_diploid is True


# Tests for is_diploid_reference()


def test_is_diploid_reference_true(diploid_fasta):
    """Test diploid detection returns True for 2-sequence file."""
    assert is_diploid_reference(diploid_fasta) is True


def test_is_diploid_reference_false_haploid(haploid_fasta):
    """Test diploid detection returns False for 1-sequence file."""
    assert is_diploid_reference(haploid_fasta) is False


def test_is_diploid_reference_false_triploid(triploid_fasta):
    """Test diploid detection returns False for 3-sequence file."""
    assert is_diploid_reference(triploid_fasta) is False


def test_is_diploid_reference_accepts_string_path(diploid_fasta):
    """Test that function accepts string path."""
    assert is_diploid_reference(str(diploid_fasta)) is True


# Tests for extract_haplotypes()


def test_extract_haplotypes_creates_two_files(diploid_fasta, temp_dir):
    """Test that haplotype extraction creates two separate files."""
    output_dir = temp_dir / "output"
    hap1_path, hap2_path = extract_haplotypes(diploid_fasta, output_dir, "test")

    assert Path(hap1_path).exists()
    assert Path(hap2_path).exists()
    assert Path(hap1_path).name == "test_hap1.fa"
    assert Path(hap2_path).name == "test_hap2.fa"


def test_extract_haplotypes_preserves_sequences(diploid_fasta, temp_dir):
    """Test that extracted haplotypes have correct sequences."""
    output_dir = temp_dir / "output"
    hap1_path, hap2_path = extract_haplotypes(diploid_fasta, output_dir, "test")

    # Read original
    original_seqs = list(SeqIO.parse(diploid_fasta, "fasta"))

    # Read extracted
    hap1_seq = next(iter(SeqIO.parse(hap1_path, "fasta")))
    hap2_seq = next(iter(SeqIO.parse(hap2_path, "fasta")))

    assert str(hap1_seq.seq) == str(original_seqs[0].seq)
    assert str(hap2_seq.seq) == str(original_seqs[1].seq)
    assert hap1_seq.id == original_seqs[0].id
    assert hap2_seq.id == original_seqs[1].id


def test_extract_haplotypes_unequal_lengths(diploid_unequal_fasta, temp_dir):
    """Test haplotype extraction with unequal lengths."""
    output_dir = temp_dir / "output"
    hap1_path, hap2_path = extract_haplotypes(diploid_unequal_fasta, output_dir, "test")

    hap1_seq = next(iter(SeqIO.parse(hap1_path, "fasta")))
    hap2_seq = next(iter(SeqIO.parse(hap2_path, "fasta")))

    assert len(hap1_seq.seq) == 20
    assert len(hap2_seq.seq) == 100


def test_extract_haplotypes_creates_output_dir(diploid_fasta, temp_dir):
    """Test that output directory is created if it doesn't exist."""
    output_dir = temp_dir / "new_dir" / "nested"
    hap1_path, hap2_path = extract_haplotypes(diploid_fasta, output_dir, "test")

    assert output_dir.exists()
    assert Path(hap1_path).exists()


def test_extract_haplotypes_uses_input_basename_if_no_base_name(diploid_fasta, temp_dir):
    """Test that input filename is used as base name when not specified."""
    output_dir = temp_dir / "output"
    hap1_path, hap2_path = extract_haplotypes(diploid_fasta, output_dir, base_name=None)

    assert "diploid_hap1.fa" in hap1_path
    assert "diploid_hap2.fa" in hap2_path


def test_extract_haplotypes_fails_on_haploid(haploid_fasta, temp_dir):
    """Test that extraction fails for non-diploid reference."""
    output_dir = temp_dir / "output"

    with pytest.raises(ValidationError, match="Expected diploid.*found 1"):
        extract_haplotypes(haploid_fasta, output_dir, "test")


def test_extract_haplotypes_fails_on_triploid(triploid_fasta, temp_dir):
    """Test that extraction fails for triploid reference."""
    output_dir = temp_dir / "output"

    with pytest.raises(ValidationError, match="Expected diploid.*found 3"):
        extract_haplotypes(triploid_fasta, output_dir, "test")


def test_extract_haplotypes_accepts_string_paths(diploid_fasta, temp_dir):
    """Test that function accepts string paths (not just Path objects)."""
    output_dir = temp_dir / "output"
    hap1_path, hap2_path = extract_haplotypes(str(diploid_fasta), str(output_dir), "test")

    assert Path(hap1_path).exists()
    assert Path(hap2_path).exists()


# Tests for validate_reference_compatibility()


def test_validate_reference_compatibility_haploid_min_1(haploid_fasta):
    """Test validation passes for haploid with min=1."""
    validate_reference_compatibility(haploid_fasta, min_sequences=1)
    # No exception = pass


def test_validate_reference_compatibility_diploid_min_2(diploid_fasta):
    """Test validation passes for diploid with min=2."""
    validate_reference_compatibility(diploid_fasta, min_sequences=2)
    # No exception = pass


def test_validate_reference_compatibility_diploid_exact_2(diploid_fasta):
    """Test validation passes for diploid with min=2, max=2."""
    validate_reference_compatibility(diploid_fasta, min_sequences=2, max_sequences=2)
    # No exception = pass


def test_validate_reference_compatibility_fails_min(haploid_fasta):
    """Test validation fails when below minimum."""
    with pytest.raises(ValidationError, match="at least 2 required"):
        validate_reference_compatibility(haploid_fasta, min_sequences=2)


def test_validate_reference_compatibility_fails_max(triploid_fasta):
    """Test validation fails when above maximum."""
    with pytest.raises(ValidationError, match="maximum 2 allowed"):
        validate_reference_compatibility(triploid_fasta, max_sequences=2)


def test_validate_reference_compatibility_accepts_string_path(diploid_fasta):
    """Test that function accepts string path."""
    validate_reference_compatibility(str(diploid_fasta), min_sequences=2, max_sequences=2)
    # No exception = pass


@pytest.mark.parametrize(
    "num_sequences,min_seq,max_seq,should_pass",
    [
        (1, 1, 1, True),  # Exact match
        (2, 1, 3, True),  # Within range
        (2, 2, 2, True),  # Exact diploid
        (1, 2, None, False),  # Below minimum
        (3, 1, 2, False),  # Above maximum
        (2, 3, None, False),  # Below minimum
    ],
)
def test_validate_reference_compatibility_parametrized(
    temp_dir, num_sequences, min_seq, max_seq, should_pass
):
    """Parametrized test for various validation scenarios."""
    # Create FASTA with specified number of sequences
    fasta_path = temp_dir / "test.fa"
    records = [
        SeqRecord(Seq("ATCG" * (i + 1)), id=f"seq{i}", description="") for i in range(num_sequences)
    ]
    SeqIO.write(records, fasta_path, "fasta")

    if should_pass:
        validate_reference_compatibility(fasta_path, min_sequences=min_seq, max_sequences=max_seq)
    else:
        with pytest.raises(ValidationError):
            validate_reference_compatibility(
                fasta_path, min_sequences=min_seq, max_sequences=max_seq
            )

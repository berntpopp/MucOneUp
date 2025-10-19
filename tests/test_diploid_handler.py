"""Comprehensive tests for diploid_handler module.

Tests cover:
- Coverage correction calculations
- Diploid simulation preparation
- Split-simulation workflow with mocked simulation
- Edge cases and error handling
"""

import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from muc_one_up.read_simulator.utils.diploid_handler import (
    DiploidSimulationResult,
    calculate_corrected_coverage,
    prepare_diploid_simulation,
    run_split_simulation,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def diploid_fasta(temp_dir):
    """Create a diploid FASTA file."""
    fasta_path = temp_dir / "diploid.fa"
    records = [
        SeqRecord(Seq("ATCG" * 50), id="haplotype_1", description=""),  # 200 bp
        SeqRecord(Seq("GCTA" * 75), id="haplotype_2", description=""),  # 300 bp
    ]
    SeqIO.write(records, fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def haploid_fasta(temp_dir):
    """Create a haploid FASTA file."""
    fasta_path = temp_dir / "haploid.fa"
    record = SeqRecord(Seq("ATCG" * 100), id="single", description="")
    SeqIO.write([record], fasta_path, "fasta")
    return fasta_path


def create_mock_fastq(path: Path, num_reads: int = 10) -> Path:
    """Helper to create a mock FASTQ file."""
    content = []
    for i in range(num_reads):
        content.append(f"@read_{i}")
        content.append("ATCGATCGATCG")
        content.append("+")
        content.append("IIIIIIIIIIII")
    path.write_text("\n".join(content) + "\n")
    return path


# Tests for calculate_corrected_coverage()


def test_calculate_corrected_coverage_diploid_default():
    """Test coverage correction for diploid split-simulation with default factor."""
    corrected = calculate_corrected_coverage(200, 0.325, is_split_simulation=True)
    # Expected: (200 / 2) / 0.325 = 307.69...
    assert corrected == pytest.approx(307.69, rel=0.01)


def test_calculate_corrected_coverage_haploid():
    """Test coverage correction for haploid (no split)."""
    corrected = calculate_corrected_coverage(200, 0.325, is_split_simulation=False)
    # Expected: 200 / 0.325 = 615.38...
    assert corrected == pytest.approx(615.38, rel=0.01)


def test_calculate_corrected_coverage_different_factors():
    """Test coverage correction with different correction factors."""
    # Factor 0.5 (50% efficiency)
    corrected1 = calculate_corrected_coverage(100, 0.5, is_split_simulation=True)
    assert corrected1 == pytest.approx(100.0, rel=0.01)  # (100/2)/0.5 = 100

    # Factor 0.25 (25% efficiency)
    corrected2 = calculate_corrected_coverage(100, 0.25, is_split_simulation=True)
    assert corrected2 == pytest.approx(200.0, rel=0.01)  # (100/2)/0.25 = 200


def test_calculate_corrected_coverage_invalid_factor_zero():
    """Test that zero correction factor raises ValueError."""
    with pytest.raises(ValueError, match="between 0 and 1"):
        calculate_corrected_coverage(200, 0.0)


def test_calculate_corrected_coverage_invalid_factor_negative():
    """Test that negative correction factor raises ValueError."""
    with pytest.raises(ValueError, match="between 0 and 1"):
        calculate_corrected_coverage(200, -0.5)


def test_calculate_corrected_coverage_invalid_factor_too_large():
    """Test that correction factor > 1 raises ValueError."""
    with pytest.raises(ValueError, match="between 0 and 1"):
        calculate_corrected_coverage(200, 1.5)


def test_calculate_corrected_coverage_invalid_coverage_zero():
    """Test that zero coverage raises ValueError."""
    with pytest.raises(ValueError, match="must be positive"):
        calculate_corrected_coverage(0, 0.325)


def test_calculate_corrected_coverage_invalid_coverage_negative():
    """Test that negative coverage raises ValueError."""
    with pytest.raises(ValueError, match="must be positive"):
        calculate_corrected_coverage(-100, 0.325)


@pytest.mark.parametrize(
    "desired,factor,is_split,expected",
    [
        (200, 0.325, True, 307.69),  # Diploid standard
        (200, 0.325, False, 615.38),  # Haploid standard
        (100, 0.5, True, 100.0),  # Diploid 50% efficiency
        (50, 0.25, True, 100.0),  # Diploid 25% efficiency
        (300, 0.3, True, 500.0),  # Diploid 30% efficiency
    ],
)
def test_calculate_corrected_coverage_parametrized(desired, factor, is_split, expected):
    """Parametrized test for various coverage calculations."""
    corrected = calculate_corrected_coverage(desired, factor, is_split_simulation=is_split)
    assert corrected == pytest.approx(expected, rel=0.01)


# Tests for prepare_diploid_simulation()


def test_prepare_diploid_simulation_creates_haplotype_files(diploid_fasta, temp_dir):
    """Test that diploid preparation creates two haplotype files."""
    output_dir = temp_dir / "sim_prep"
    hap1, hap2, prep_dir = prepare_diploid_simulation(diploid_fasta, output_dir)

    assert Path(hap1).exists()
    assert Path(hap2).exists()
    assert "hap1.fa" in hap1
    assert "hap2.fa" in hap2
    assert prep_dir == str(output_dir)


def test_prepare_diploid_simulation_preserves_sequences(diploid_fasta, temp_dir):
    """Test that haplotype sequences match original diploid."""
    output_dir = temp_dir / "sim_prep"
    hap1, hap2, _ = prepare_diploid_simulation(diploid_fasta, output_dir)

    # Read original
    original = list(SeqIO.parse(diploid_fasta, "fasta"))

    # Read extracted
    hap1_seq = next(iter(SeqIO.parse(hap1, "fasta")))
    hap2_seq = next(iter(SeqIO.parse(hap2, "fasta")))

    assert str(hap1_seq.seq) == str(original[0].seq)
    assert str(hap2_seq.seq) == str(original[1].seq)


def test_prepare_diploid_simulation_custom_base_name(diploid_fasta, temp_dir):
    """Test using custom base name for output files."""
    output_dir = temp_dir / "sim_prep"
    hap1, hap2, _ = prepare_diploid_simulation(diploid_fasta, output_dir, base_name="custom")

    assert "custom_hap1.fa" in hap1
    assert "custom_hap2.fa" in hap2


def test_prepare_diploid_simulation_fails_on_haploid(haploid_fasta, temp_dir):
    """Test that preparation fails for non-diploid reference."""
    output_dir = temp_dir / "sim_prep"

    with pytest.raises(ValueError, match="Expected diploid.*found 1"):
        prepare_diploid_simulation(haploid_fasta, output_dir)


def test_prepare_diploid_simulation_accepts_string_paths(diploid_fasta, temp_dir):
    """Test that function accepts string paths."""
    output_dir = temp_dir / "sim_prep"
    hap1, hap2, prep_dir = prepare_diploid_simulation(str(diploid_fasta), str(output_dir))

    assert Path(hap1).exists()
    assert Path(hap2).exists()


# Tests for run_split_simulation()


def test_run_split_simulation_basic_workflow(diploid_fasta, temp_dir):
    """Test complete split-simulation workflow with mocked simulation function."""
    output_fastq = temp_dir / "merged.fastq"

    # Create mock simulation function
    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        # Create a mock FASTQ based on which haplotype
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        num_reads = 10 if "hap1" in output_prefix else 15
        create_mock_fastq(fastq_path, num_reads)
        return str(fastq_path)

    # Run split simulation
    result = run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200, "threads": 4},
        output_fastq,
        correction_factor=0.325,
        seed=42,
    )

    # Verify result
    assert isinstance(result, DiploidSimulationResult)
    assert Path(result.merged_fastq).exists()
    assert result.reads_hap1 == 10
    assert result.reads_hap2 == 15


def test_run_split_simulation_uses_corrected_coverage(diploid_fasta, temp_dir):
    """Test that split simulation passes corrected coverage to simulation function."""
    output_fastq = temp_dir / "merged.fastq"
    coverage_used = []

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        coverage_used.append(coverage)
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200, "threads": 4},
        output_fastq,
        correction_factor=0.325,
    )

    # Both haplotypes should get same corrected coverage
    assert len(coverage_used) == 2
    expected_corrected = calculate_corrected_coverage(200, 0.325, is_split_simulation=True)
    assert coverage_used[0] == pytest.approx(expected_corrected, rel=0.01)
    assert coverage_used[1] == pytest.approx(expected_corrected, rel=0.01)


def test_run_split_simulation_uses_different_seeds(diploid_fasta, temp_dir):
    """Test that haplotypes get different seeds (seed and seed+1)."""
    output_fastq = temp_dir / "merged.fastq"
    seeds_used = []

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        seeds_used.append(seed)
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200},
        output_fastq,
        seed=42,
    )

    assert len(seeds_used) == 2
    assert seeds_used[0] == 42
    assert seeds_used[1] == 43


def test_run_split_simulation_no_seed_passes_none(diploid_fasta, temp_dir):
    """Test that no seed results in None being passed."""
    output_fastq = temp_dir / "merged.fastq"
    seeds_used = []

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        seeds_used.append(seed)
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200},
        output_fastq,
        seed=None,
    )

    assert len(seeds_used) == 2
    assert seeds_used[0] is None
    assert seeds_used[1] is None


def test_run_split_simulation_passes_extra_params(diploid_fasta, temp_dir):
    """Test that extra parameters are passed through to simulation function."""
    output_fastq = temp_dir / "merged.fastq"
    received_params = []

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        received_params.append(params)
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200, "threads": 8, "custom_param": "value"},
        output_fastq,
    )

    # Both calls should receive extra params
    assert len(received_params) == 2
    assert received_params[0]["threads"] == 8
    assert received_params[0]["custom_param"] == "value"
    assert received_params[1]["threads"] == 8
    assert received_params[1]["custom_param"] == "value"


def test_run_split_simulation_merged_fastq_has_all_reads(diploid_fasta, temp_dir):
    """Test that merged FASTQ contains reads from both haplotypes."""
    output_fastq = temp_dir / "merged.fastq"

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        num_reads = 10 if "hap1" in output_prefix else 15
        create_mock_fastq(fastq_path, num_reads)
        return str(fastq_path)

    result = run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200},
        output_fastq,
    )

    # Check merged file has correct total
    from muc_one_up.read_simulator.utils.fastq_utils import count_fastq_reads

    total_reads = count_fastq_reads(result.merged_fastq)
    assert total_reads == 25  # 10 + 15


def test_run_split_simulation_creates_output_directory(diploid_fasta, temp_dir):
    """Test that output directory is created if it doesn't exist."""
    output_fastq = temp_dir / "nested" / "dir" / "merged.fastq"

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    result = run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200},
        output_fastq,
    )

    assert Path(result.merged_fastq).exists()
    assert Path(result.merged_fastq).parent.exists()


def test_run_split_simulation_accepts_string_paths(diploid_fasta, temp_dir):
    """Test that function accepts string paths."""
    output_fastq = temp_dir / "merged.fastq"

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    result = run_split_simulation(
        str(diploid_fasta),
        mock_sim_func,
        {"coverage": 200},
        str(output_fastq),
    )

    assert Path(result.merged_fastq).exists()


def test_run_split_simulation_result_structure(diploid_fasta, temp_dir):
    """Test that result has correct structure and data types."""
    output_fastq = temp_dir / "merged.fastq"

    def mock_sim_func(reference_fasta, output_prefix, coverage, seed, **params):
        fastq_path = Path(f"{output_prefix}_aligned_reads.fastq")
        create_mock_fastq(fastq_path, 10)
        return str(fastq_path)

    result = run_split_simulation(
        diploid_fasta,
        mock_sim_func,
        {"coverage": 200},
        output_fastq,
    )

    # Check all fields are present and correct types
    assert isinstance(result.merged_fastq, str)
    assert isinstance(result.hap1_fastq, str)
    assert isinstance(result.hap2_fastq, str)
    assert isinstance(result.hap1_reference, str)
    assert isinstance(result.hap2_reference, str)
    assert isinstance(result.reads_hap1, int)
    assert isinstance(result.reads_hap2, int)

    # Check merged FASTQ exists (others are in temp dir and cleaned up)
    assert Path(result.merged_fastq).exists()

    # Check reference paths have expected format (temp files are cleaned up)
    assert "hap1.fa" in result.hap1_reference
    assert "hap2.fa" in result.hap2_reference

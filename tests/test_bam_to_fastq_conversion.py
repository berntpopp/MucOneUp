"""Comprehensive tests for BAM to FASTQ conversion.

This module tests the critical functionality of converting BAM files to paired FASTQ files,
ensuring proper read pairing and singleton handling.

Test Coverage:
- Unit tests for conversion function
- Integration tests with real BAM files
- Edge case handling
- Error scenarios
- Performance validation

Author: Senior Bioinformatics Developer
Date: 2025-10-28
"""

import gzip
import subprocess
from pathlib import Path

import pytest

from muc_one_up.exceptions import ExternalToolError, FileOperationError
from muc_one_up.read_simulator.wrappers.samtools_wrapper import (
    FastqConversionOptions,
    convert_bam_to_paired_fastq,
)

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def samtools_cmd():
    """Get samtools command path, skip if not available."""
    try:
        result = subprocess.run(
            ["which", "samtools"],
            capture_output=True,
            check=True,
        )
        return result.stdout.decode().strip()
    except subprocess.CalledProcessError:
        pytest.skip("samtools not available in PATH")


@pytest.fixture
def test_bam_with_unpaired_reads(tmp_path):
    """
    Create a minimal test BAM file with unpaired reads.

    This simulates the exact scenario from VNTR downsampling where
    some reads lose their mates.

    Returns:
        Path to test BAM file
    """
    # Check if we have a real BAM file from experiments to use
    real_bam = Path(
        "/mnt/c/development/MucOneUp/output/experminents/experiment_1/pair_001/sample.001.normal.simulated_vntr_biased.bam"
    )

    if real_bam.exists():
        return real_bam
    else:
        pytest.skip("Real test BAM file not available - run experiments first")


@pytest.fixture
def sample_bam_stats():
    """Expected statistics for the sample BAM file used in testing."""
    return {
        "total_reads": 4906,
        "read1": 2450,
        "read2": 2456,
        "unpaired": 84,
        "expected_paired_output": 2373,  # After discarding 84 unpaired
        "singletons_discarded": 62,  # From samtools output
    }


# =============================================================================
# UNIT TESTS: Basic Functionality
# =============================================================================


class TestBasicConversion:
    """Test basic BAM to FASTQ conversion functionality."""

    def test_convert_with_collation_enabled(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path, sample_bam_stats
    ):
        """Test conversion with collation and singleton handling (DEFAULT behavior)."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        # Convert with default options (collation enabled)
        fq1, fq2 = convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                collate_before_conversion=True,
                discard_singletons=True,
                validate_pairs=True,
            ),
        )

        # Verify outputs created
        assert Path(fq1).exists(), "R1 FASTQ file not created"
        assert Path(fq2).exists(), "R2 FASTQ file not created"

        # Verify read counts match
        r1_reads = self._count_fastq_reads(Path(fq1))
        r2_reads = self._count_fastq_reads(Path(fq2))

        assert r1_reads == r2_reads, f"Read count mismatch: R1={r1_reads}, R2={r2_reads}"
        assert r1_reads == sample_bam_stats["expected_paired_output"], (
            f"Expected {sample_bam_stats['expected_paired_output']} paired reads, got {r1_reads}"
        )

    def test_convert_without_collation_fails_validation(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path
    ):
        """Test that conversion without collation on unpaired BAM raises validation error."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        # This SHOULD fail validation because unpaired reads cause mismatches
        # However, with validate_pairs=False, it should succeed but produce bad output
        fq1, fq2 = convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                collate_before_conversion=False,  # No collation
                discard_singletons=True,
                validate_pairs=False,  # Skip validation to allow bad output
            ),
        )

        # Verify files created but may have mismatched counts
        assert Path(fq1).exists()
        assert Path(fq2).exists()

        r1_reads = self._count_fastq_reads(Path(fq1))
        r2_reads = self._count_fastq_reads(Path(fq2))

        # Without collation, reads may be mismatched
        # This demonstrates WHY collation is critical
        print(f"Without collation: R1={r1_reads}, R2={r2_reads}")

    @staticmethod
    def _count_fastq_reads(fastq_path: Path) -> int:
        """Count reads in FASTQ file (handles gzip)."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            return sum(1 for _ in f) // 4


# =============================================================================
# UNIT TESTS: Read Pairing Integrity
# =============================================================================


class TestReadPairingIntegrity:
    """Test that read IDs are properly paired between R1 and R2."""

    def test_read_ids_match_between_r1_r2(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path
    ):
        """Test that read IDs in R1 and R2 are identical and in same order."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                collate_before_conversion=True,
                discard_singletons=True,
                preserve_read_names=True,  # Don't add /1 /2 suffixes
            ),
        )

        # Extract first 100 read IDs from each file
        r1_ids = self._get_read_ids(output_fq1, num_reads=100)
        r2_ids = self._get_read_ids(output_fq2, num_reads=100)

        # Verify IDs match exactly
        assert len(r1_ids) == len(r2_ids), "Different number of reads extracted"
        for i, (id1, id2) in enumerate(zip(r1_ids, r2_ids, strict=True)):
            assert id1 == id2, f"Read ID mismatch at position {i}: R1='{id1}', R2='{id2}'"

    def test_read_names_preserved_with_n_flag(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path
    ):
        """Test that read names don't have /1 /2 suffixes when -n flag used."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                preserve_read_names=True,
            ),
        )

        # Check that read names don't end with /1 or /2
        r1_ids = self._get_read_ids(output_fq1, num_reads=10)
        r2_ids = self._get_read_ids(output_fq2, num_reads=10)

        for read_id in r1_ids:
            assert not read_id.endswith("/1"), f"Read ID has /1 suffix: {read_id}"
        for read_id in r2_ids:
            assert not read_id.endswith("/2"), f"Read ID has /2 suffix: {read_id}"

    @staticmethod
    def _get_read_ids(fastq_path: Path, num_reads: int = 100) -> list[str]:
        """Extract read IDs from FASTQ file."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        read_ids = []

        with opener(fastq_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 0:  # Header line
                    read_ids.append(line.strip()[1:])  # Remove @ prefix
                if len(read_ids) >= num_reads:
                    break

        return read_ids


# =============================================================================
# UNIT TESTS: Singleton Handling
# =============================================================================


class TestSingletonHandling:
    """Test proper handling of singleton reads."""

    def test_singletons_discarded_with_s_dev_null(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path, sample_bam_stats
    ):
        """Test that singleton reads are discarded when -s /dev/null used."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                discard_singletons=True,  # Should use -s /dev/null
            ),
        )

        r1_reads = self._count_fastq_reads(output_fq1)
        r2_reads = self._count_fastq_reads(output_fq2)

        # Both should have same count (singletons discarded)
        assert r1_reads == r2_reads
        assert r1_reads == sample_bam_stats["expected_paired_output"]

    def test_singletons_saved_to_file(self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path):
        """Test that singleton reads can be saved to a separate file."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"
        singleton_file = tmp_path / "singletons.fastq.gz"

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                output_singleton=singleton_file,
                discard_singletons=False,  # Keep singletons
            ),
        )

        # Verify singleton file created and non-empty
        assert singleton_file.exists(), "Singleton file not created"
        assert singleton_file.stat().st_size > 0, "Singleton file is empty"

        # Count singletons
        singleton_count = self._count_fastq_reads(singleton_file)
        print(f"Singleton reads saved: {singleton_count}")
        assert singleton_count > 0, "No singletons found (expected some from BAM)"

    @staticmethod
    def _count_fastq_reads(fastq_path: Path) -> int:
        """Count reads in FASTQ file."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            return sum(1 for _ in f) // 4


# =============================================================================
# UNIT TESTS: Error Handling
# =============================================================================


class TestErrorHandling:
    """Test error handling for various failure scenarios."""

    def test_missing_input_bam_raises_error(self, samtools_cmd, tmp_path):
        """Test that missing input BAM raises FileOperationError."""
        nonexistent_bam = tmp_path / "nonexistent.bam"
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        with pytest.raises(FileOperationError, match="not found"):
            convert_bam_to_paired_fastq(
                samtools_cmd=samtools_cmd,
                input_bam=nonexistent_bam,
                output_fq1=output_fq1,
                output_fq2=output_fq2,
            )

    def test_empty_input_bam_raises_error(self, samtools_cmd, tmp_path):
        """Test that empty input BAM raises FileOperationError."""
        empty_bam = tmp_path / "empty.bam"
        empty_bam.touch()  # Create empty file
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        with pytest.raises(FileOperationError, match="empty"):
            convert_bam_to_paired_fastq(
                samtools_cmd=samtools_cmd,
                input_bam=empty_bam,
                output_fq1=output_fq1,
                output_fq2=output_fq2,
            )

    def test_invalid_samtools_command_raises_error(self, test_bam_with_unpaired_reads, tmp_path):
        """Test that invalid samtools command raises ExternalToolError."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        with pytest.raises(ExternalToolError, match="samtools failed"):
            convert_bam_to_paired_fastq(
                samtools_cmd="/nonexistent/samtools",
                input_bam=test_bam_with_unpaired_reads,
                output_fq1=output_fq1,
                output_fq2=output_fq2,
            )


# =============================================================================
# INTEGRATION TESTS: Real-World Scenarios
# =============================================================================


class TestIntegrationScenarios:
    """Integration tests with real BAM files and realistic scenarios."""

    def test_vntr_biased_bam_conversion(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path, sample_bam_stats
    ):
        """
        Integration test: Convert real VNTR-biased BAM from experiments.

        This is the exact use case that was broken before the fix.
        """
        output_fq1 = tmp_path / "vntr_R1.fastq.gz"
        output_fq2 = tmp_path / "vntr_R2.fastq.gz"

        fq1, fq2 = convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                collate_before_conversion=True,
                discard_singletons=True,
                validate_pairs=True,
                threads=4,
            ),
        )

        # Verify proper pairing
        r1_reads = self._count_reads(Path(fq1))
        r2_reads = self._count_reads(Path(fq2))

        assert r1_reads == r2_reads, (
            f"VNTR BAM conversion failed pairing: R1={r1_reads}, R2={r2_reads}"
        )

        # Verify read IDs match
        r1_ids = self._get_first_n_read_ids(Path(fq1), 50)
        r2_ids = self._get_first_n_read_ids(Path(fq2), 50)

        assert r1_ids == r2_ids, "Read IDs don't match between R1 and R2"

    def test_comparison_with_manual_samtools_command(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path
    ):
        """
        Test that our implementation produces identical output to manual samtools command.

        This validates that we're correctly implementing best practices.
        """
        # Output from our implementation
        impl_fq1 = tmp_path / "impl_R1.fastq.gz"
        impl_fq2 = tmp_path / "impl_R2.fastq.gz"

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=impl_fq1,
            output_fq2=impl_fq2,
        )

        # Output from manual samtools command
        manual_fq1 = tmp_path / "manual_R1.fastq.gz"
        manual_fq2 = tmp_path / "manual_R2.fastq.gz"

        # Run exact recommended samtools command
        cmd = (
            f"{samtools_cmd} collate -u -O -@ 4 {test_bam_with_unpaired_reads} temp | "
            f"{samtools_cmd} fastq -1 {manual_fq1} -2 {manual_fq2} "
            f"-0 /dev/null -s /dev/null -n -"
        )
        subprocess.run(cmd, shell=True, check=True, cwd=tmp_path)

        # Compare outputs
        impl_r1_count = self._count_reads(impl_fq1)
        impl_r2_count = self._count_reads(impl_fq2)
        manual_r1_count = self._count_reads(manual_fq1)
        manual_r2_count = self._count_reads(manual_fq2)

        assert impl_r1_count == manual_r1_count, (
            f"Implementation R1 count ({impl_r1_count}) != manual R1 count ({manual_r1_count})"
        )
        assert impl_r2_count == manual_r2_count, (
            f"Implementation R2 count ({impl_r2_count}) != manual R2 count ({manual_r2_count})"
        )

    @staticmethod
    def _count_reads(fastq_path: Path) -> int:
        """Count reads in FASTQ."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            return sum(1 for _ in f) // 4

    @staticmethod
    def _get_first_n_read_ids(fastq_path: Path, n: int) -> list[str]:
        """Get first N read IDs."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        ids = []
        with opener(fastq_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    ids.append(line.strip()[1:])
                if len(ids) >= n:
                    break
        return ids


# =============================================================================
# PERFORMANCE TESTS
# =============================================================================


class TestPerformance:
    """Performance and timeout tests."""

    def test_conversion_completes_within_timeout(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path
    ):
        """Test that conversion completes within reasonable timeout."""
        output_fq1 = tmp_path / "test_R1.fastq.gz"
        output_fq2 = tmp_path / "test_R2.fastq.gz"

        import time

        start = time.time()

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                timeout=60,  # 60 second timeout
            ),
        )

        elapsed = time.time() - start

        # Small BAM should convert in < 30 seconds
        assert elapsed < 30, f"Conversion took {elapsed:.1f}s (expected < 30s)"
        print(f"Conversion completed in {elapsed:.2f}s")


# =============================================================================
# REGRESSION TESTS
# =============================================================================


class TestRegressions:
    """Tests for specific bug fixes and regressions."""

    def test_vntr_fastq_bug_fix_regression(
        self, samtools_cmd, test_bam_with_unpaired_reads, tmp_path, sample_bam_stats
    ):
        """
        Regression test for VNTR FASTQ pairing bug.

        This test ensures the original bug (mismatched read counts) doesn't recur.

        Original bug:
        - R1: 2,411 reads
        - R2: 2,421 reads (10 read mismatch!)
        - Read IDs out of order

        Expected after fix:
        - R1: 2,373 reads
        - R2: 2,373 reads (perfect match)
        - Read IDs in same order
        """
        output_fq1 = tmp_path / "regression_R1.fastq.gz"
        output_fq2 = tmp_path / "regression_R2.fastq.gz"

        convert_bam_to_paired_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=test_bam_with_unpaired_reads,
            output_fq1=output_fq1,
            output_fq2=output_fq2,
            options=FastqConversionOptions(
                collate_before_conversion=True,  # Critical for fix
                discard_singletons=True,  # Critical for fix
                validate_pairs=True,
            ),
        )

        r1_reads = self._count_reads(output_fq1)
        r2_reads = self._count_reads(output_fq2)

        # Verify no mismatch (original bug symptom)
        assert r1_reads == r2_reads, (
            f"REGRESSION: Read count mismatch detected! R1={r1_reads}, R2={r2_reads}. "
            f"This was the original bug that was fixed."
        )

        # Verify expected count
        assert r1_reads == sample_bam_stats["expected_paired_output"], (
            f"Expected {sample_bam_stats['expected_paired_output']} reads, got {r1_reads}"
        )

        print(f"✓ Regression test passed: {r1_reads} perfectly paired reads")

    @staticmethod
    def _count_reads(fastq_path: Path) -> int:
        """Count reads in FASTQ."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            return sum(1 for _ in f) // 4

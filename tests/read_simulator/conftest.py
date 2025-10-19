"""Shared fixtures for read simulator tests.

This module provides common fixtures for testing read simulator wrappers
and pipelines. Following best practices, we mock only at system boundaries
(subprocess calls) to test OUR code's logic.
"""

from unittest.mock import Mock

import pytest


@pytest.fixture
def mock_subprocess_popen(mocker):
    """Mock subprocess.Popen for testing process chaining (BWA | samtools).

    Returns a factory function that creates mock process objects with
    configurable return codes and output.

    Usage:
        processes = mock_subprocess_popen(return_codes=[0, 0])
        # First call returns processes[0], second call returns processes[1]
    """

    def _create_mocks(return_codes=None, stdouts=None, stderrs=None):
        """Create mock Popen processes with specified return codes and outputs."""
        if return_codes is None:
            return_codes = [0]
        if stdouts is None:
            stdouts = [b""] * len(return_codes)
        if stderrs is None:
            stderrs = [b""] * len(return_codes)

        mock_processes = []
        for returncode, stdout, stderr in zip(return_codes, stdouts, stderrs, strict=False):
            mock_process = Mock()
            mock_process.returncode = returncode
            mock_process.stdout = Mock()
            mock_process.stdout.close = Mock()
            mock_process.communicate = Mock(return_value=(stdout, stderr))

            mock_processes.append(mock_process)

        # Patch subprocess.Popen to return mocks in sequence
        mock_popen = mocker.patch("subprocess.Popen", side_effect=mock_processes)

        return {"popen": mock_popen, "processes": mock_processes}

    return _create_mocks


@pytest.fixture
def mock_subprocess_run(mocker):
    """Mock subprocess.run for testing simple command execution.

    Returns a mock with configurable return code, stdout, and stderr.

    Usage:
        result = mock_subprocess_run(returncode=0, stderr=b"")
        # Calls to subprocess.run will return this mock
    """

    def _create_mock(returncode=0, stdout=b"", stderr=b""):
        """Create mock CompletedProcess with specified values."""
        mock_result = Mock()
        mock_result.returncode = returncode
        mock_result.stdout = stdout
        mock_result.stderr = stderr

        mock_run = mocker.patch("subprocess.run", return_value=mock_result)

        return {"run": mock_run, "result": mock_result}

    return _create_mock


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file for testing.

    Returns:
        Path to a temporary FASTA file with one sequence
    """
    fasta = tmp_path / "sample.fa"
    fasta.write_text(">chr1\nACGTACGTACGTACGT\n")
    return fasta


@pytest.fixture
def sample_fastq_files(tmp_path):
    """Create sample paired FASTQ files for testing.

    Returns:
        Tuple of (read1_path, read2_path)
    """
    r1 = tmp_path / "read1.fq"
    r2 = tmp_path / "read2.fq"

    r1.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n")
    r2.write_text("@read1\nTGCATGCA\n+\nIIIIIIII\n")

    return (r1, r2)


@pytest.fixture
def sample_bam_file(tmp_path):
    """Create a fake BAM file for testing.

    Returns:
        Path to a temporary file with BAM-like content
    """
    bam = tmp_path / "sample.bam"
    # Write minimal BAM header bytes
    bam.write_bytes(b"BAM\x01FAKE_BAM_DATA_FOR_TESTING")
    return bam


@pytest.fixture
def tools_dict():
    """Standard tools dictionary for testing.

    Returns:
        Dict mapping tool names to command strings
    """
    return {
        "bwa": "bwa",
        "samtools": "samtools",
        "reseq": "reseq",
        "nanosim": "nanosim",
        "minimap2": "minimap2",
        "pblat": "pblat",
        "faToTwoBit": "faToTwoBit",
    }

"""Shared pytest fixtures for MucOneUp tests.

This module provides reusable fixtures for testing MucOneUp functionality,
following DRY principles to avoid test code duplication.
"""

import json
from pathlib import Path
from typing import Any

import pytest

# ============================================================================
# Configuration Fixtures
# ============================================================================


@pytest.fixture
def minimal_config() -> dict[str, Any]:
    """Minimal valid configuration for basic tests.

    This fixture provides a complete but minimal configuration that satisfies
    all schema requirements. Use this as a base for most tests.
    """
    return {
        "repeats": {
            "1": "ATCGATCGATCG",
            "2": "GCTAGCTAGCTA",
            "X": "AAACCCGGGTTT",
            "A": "TTTTGGGGCCCC",
            "B": "AAAATTTTGGGG",
            "C": "CCCCCCCCCCCC",
            "6": "GGGCCCAAATTT",
            "6p": "TTTAAACCCGGG",
            "7": "AAATAAA",
            "8": "TTTCTTT",
            "9": "GGGAGGG",
        },
        "constants": {
            "hg38": {
                "left": "AAAAAAAAAA" * 10,  # 100 bp
                "right": "TTTTTTTTTT" * 10,  # 100 bp
                "vntr_start": 100,
                "vntr_end": 500,
            }
        },
        "probabilities": {
            "1": {"2": 1.0},
            "2": {"X": 0.5, "A": 0.5},
            "X": {"B": 0.5, "6": 0.5},
            "A": {"B": 0.5, "6": 0.5},
            "B": {"6": 0.5, "6p": 0.5},
            "C": {"X": 1.0},
            "6": {"7": 1.0},
            "6p": {"7": 1.0},
            "7": {"8": 1.0},
            "8": {"9": 1.0},
            "9": {},
        },
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
            "median_repeats": 30,
        },
        "mutations": {
            "dupC": {
                "allowed_repeats": ["X", "A", "B"],
                "strict_mode": False,
                "changes": [
                    {
                        "type": "insert",
                        "position": 1,
                        "sequence": "CCCCCCCC",
                    }
                ],
            },
            "delA": {
                "allowed_repeats": ["A"],
                "strict_mode": True,
                "changes": [
                    {
                        "type": "delete",
                        "position": 0,
                        "count": 5,
                    }
                ],
            },
        },
        "tools": {
            "samtools": "samtools",
            "bwa": "bwa",
            "reseq": "reseq",
            "minimap2": "minimap2",
        },
        "read_simulation": {
            "simulator": "illumina",
            "human_reference": "/path/to/ref.fa",
            "threads": 2,
            "fragment_size_mean": 300,
            "fragment_size_sd": 50,
            "coverage": 30,
        },
        "reference_assembly": "hg38",
    }


@pytest.fixture
def temp_config_file(tmp_path: Path, minimal_config: dict[str, Any]) -> Path:
    """Create temporary config file.

    Args:
        tmp_path: pytest temporary directory
        minimal_config: Configuration dictionary

    Returns:
        Path to temporary config file
    """
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(minimal_config, indent=2))
    return config_file


# ============================================================================
# File Fixtures
# ============================================================================


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    """Create sample FASTA file with two haplotypes.

    Returns:
        Path to temporary FASTA file
    """
    fasta = tmp_path / "sample.fa"
    fasta.write_text(
        ">haplotype_1\n" "ATCGATCGATCGATCGATCGATCG\n" ">haplotype_2\n" "GCTAGCTAGCTAGCTAGCTAGCTA\n"
    )
    return fasta


@pytest.fixture
def sample_structure_file(tmp_path: Path) -> Path:
    """Create sample structure file.

    Returns:
        Path to temporary structure file
    """
    structure = tmp_path / "structure.txt"
    structure.write_text("haplotype_1\t1-2-X-B-6-7-8-9\n" "haplotype_2\t1-2-A-B-6p-7-8-9\n")
    return structure


@pytest.fixture
def sample_structure_file_with_mutations(tmp_path: Path) -> Path:
    """Create sample structure file with mutation markers.

    Returns:
        Path to temporary structure file with mutations
    """
    structure = tmp_path / "structure_mutated.txt"
    structure.write_text(
        "# Mutation Applied: dupC (Targets: [(1, 25)])\n"
        "haplotype_1\t1-2-X-Bm-6-7-8-9\n"
        "haplotype_2\t1-2-A-B-6p-7-8-9\n"
    )
    return structure


@pytest.fixture
def sample_snp_file(tmp_path: Path) -> Path:
    """Create sample SNP TSV file.

    Format: haplotype(1-based) position(0-based) ref alt

    Returns:
        Path to temporary SNP file
    """
    snp_file = tmp_path / "snps.tsv"
    snp_file.write_text(
        "1\t110\tA\tG\n"  # Haplotype 1, position 110, A→G
        "2\t120\tC\tT\n"  # Haplotype 2, position 120, C→T
        "1\t130\tG\tC\n"  # Haplotype 1, position 130, G→C
    )
    return snp_file


# ============================================================================
# Output Directory Fixtures
# ============================================================================


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    """Create temporary output directory.

    Returns:
        Path to temporary output directory
    """
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    return out_dir


# ============================================================================
# Sequence Fixtures (Bioinformatics)
# ============================================================================


@pytest.fixture
def valid_dna_sequence() -> str:
    """Valid DNA sequence for testing (240 bp).

    Returns:
        DNA sequence string containing only ACGT
    """
    return "ATCGATCGATCGATCGATCGATCG" * 10


@pytest.fixture
def invalid_dna_sequence() -> str:
    """Invalid DNA sequence containing non-ACGTN characters.

    Returns:
        DNA sequence with invalid characters (X, Y, Z)
    """
    return "ATCGXYZATCG"


@pytest.fixture
def sample_haplotypes() -> list[tuple[str, list[str]]]:
    """Sample haplotype results for testing.

    Returns:
        List of (haplotype_name, repeat_chain) tuples
    """
    return [
        ("haplotype_1", ["1", "2", "X", "B", "6", "7", "8", "9"]),
        ("haplotype_2", ["1", "2", "A", "B", "6p", "7", "8", "9"]),
    ]


@pytest.fixture
def sample_haplotype_sequences(minimal_config: dict[str, Any]) -> list[tuple[str, str]]:
    """Sample haplotype sequences (name, sequence) for testing.

    Args:
        minimal_config: Configuration fixture

    Returns:
        List of (haplotype_name, dna_sequence) tuples
    """
    # Build simple sequences from repeats
    config = minimal_config
    constants = config["constants"]["hg38"]
    repeats = config["repeats"]

    seq1 = (
        constants["left"]
        + repeats["1"]
        + repeats["2"]
        + repeats["X"]
        + repeats["B"]
        + repeats["6"]
        + repeats["7"]
        + repeats["8"]
        + repeats["9"]
        + constants["right"]
    )

    seq2 = (
        constants["left"]
        + repeats["1"]
        + repeats["2"]
        + repeats["A"]
        + repeats["B"]
        + repeats["6p"]
        + repeats["7"]
        + repeats["8"]
        + repeats["9"]
        + constants["right"]
    )

    return [("haplotype_1", seq1), ("haplotype_2", seq2)]


# ============================================================================
# Mock Fixtures
# ============================================================================


@pytest.fixture
def mock_external_tools(monkeypatch):
    """Mock external tool execution for faster tests.

    This fixture mocks subprocess.run to avoid calling actual external tools
    like samtools, bwa, etc. during unit tests.

    Args:
        monkeypatch: pytest monkeypatch fixture
    """
    from unittest.mock import Mock

    def mock_subprocess_run(cmd, *args, **kwargs):
        """Mock subprocess.run for external tools."""
        result = Mock()
        result.returncode = 0
        result.stdout = b"Mock output"
        result.stderr = b""
        return result

    import subprocess

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)


# ============================================================================
# Cleanup
# ============================================================================


@pytest.fixture(autouse=True)
def cleanup_logging():
    """Auto-cleanup of logging handlers after each test.

    This prevents log handler accumulation across tests.
    """
    import logging

    yield

    # Remove all handlers from root logger
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        handler.close()
        root_logger.removeHandler(handler)

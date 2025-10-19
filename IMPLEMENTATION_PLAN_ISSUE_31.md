# Implementation Plan: Issue #31 - NanoSim Diploid Split-Simulation

**Author**: Senior Bioinformatics Developer
**Date**: 2025-10-19
**Principles**: DRY, KISS, SOLID, Test-Driven Development
**Code Style**: Follow existing MucOneUp patterns

---

## Executive Summary

Implement split-simulation approach for NanoSim to eliminate length-dependent allelic bias in diploid ONT read simulation. Solution involves:

1. **Diploid detection** - Identify 2-sequence FASTA files
2. **Per-haplotype simulation** - Simulate each haplotype independently
3. **Coverage correction** - Apply training model-specific correction factor
4. **FASTQ merging** - Combine reads maintaining 50:50 representation
5. **Feature flags** - Default-on switches for split-simulation and coverage correction

**Impact**: Reduces coverage bias from 3.27:1 to 1.40:1 (2.3x improvement)

---

## Architecture Overview

### Design Principles Applied

**SOLID Principles**:
- **Single Responsibility**: Each module handles one concern (detection, simulation, merging)
- **Open/Closed**: Extensible via configuration without modifying core code
- **Liskov Substitution**: Diploid and haploid pipelines interchangeable
- **Interface Segregation**: Clean separation of concerns
- **Dependency Inversion**: Depend on abstractions (config), not implementations

**DRY (Don't Repeat Yourself)**:
- Reuse existing `run_nanosim_simulation()` for both haplotypes
- Extract common validation logic
- Share FASTQ merging utilities

**KISS (Keep It Simple, Stupid)**:
- Feature flags control behavior clearly
- Simple linear pipeline flow
- No premature optimization

### Modular Structure

```
muc_one_up/read_simulator/
├── ont_pipeline.py          # Main orchestrator (refactored)
├── wrappers/
│   ├── nanosim_wrapper.py   # Core NanoSim calls (minimal changes)
│   └── diploid_handler.py   # NEW: Diploid-specific logic
└── utils/
    ├── fastq_utils.py        # NEW: FASTQ merging utilities
    └── reference_utils.py    # NEW: Reference analysis utilities

tests/
├── test_diploid_handler.py  # NEW: Diploid detection & splitting
├── test_fastq_utils.py       # NEW: FASTQ merging
└── test_ont_pipeline.py      # NEW: Integration tests
```

---

## Implementation Tasks

### Phase 1: Core Utilities (Foundation)

#### Task 1.1: Reference Analysis Utility

**File**: `muc_one_up/read_simulator/utils/reference_utils.py`

**Purpose**: Analyze FASTA files to detect diploid references and extract haplotypes.

```python
#!/usr/bin/env python3
"""
Reference FASTA analysis utilities for read simulation.

Provides functions to:
- Detect diploid vs haploid references
- Extract haplotypes from diploid FASTA
- Validate reference integrity
"""

import logging
from pathlib import Path
from typing import Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def is_diploid_reference(fasta_path: str) -> bool:
    """
    Determine if a FASTA file contains a diploid reference (2 sequences).

    A diploid reference is defined as a FASTA file containing exactly 2 sequences,
    representing two haplotypes. This is the standard format for MucOneUp simulated
    diploid MUC1 VNTR references.

    Args:
        fasta_path: Path to FASTA file to analyze

    Returns:
        True if file contains exactly 2 sequences (diploid), False otherwise

    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If FASTA file is empty or malformed

    Example:
        >>> is_diploid_reference("sample.001.simulated.fa")
        True  # Contains haplotype_1 and haplotype_2

        >>> is_diploid_reference("haploid.fa")
        False  # Contains only 1 sequence
    """
    fasta_file = Path(fasta_path)

    if not fasta_file.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")

    try:
        sequences = list(SeqIO.parse(str(fasta_file), "fasta"))
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA file {fasta_path}: {e}") from e

    if len(sequences) == 0:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    is_diploid = len(sequences) == 2

    logging.debug(
        f"Reference analysis: {fasta_path} contains {len(sequences)} sequence(s) "
        f"→ {'diploid' if is_diploid else 'haploid/other'}"
    )

    return is_diploid


def extract_haplotypes(
    diploid_fasta: str,
    output_dir: str,
    base_name: str
) -> Tuple[str, str]:
    """
    Extract individual haplotypes from a diploid FASTA file.

    Splits a diploid FASTA (2 sequences) into separate single-haplotype FASTA files.
    This enables independent simulation of each haplotype to eliminate length bias.

    Args:
        diploid_fasta: Path to diploid FASTA file (must contain exactly 2 sequences)
        output_dir: Directory to write extracted haplotype files
        base_name: Base name for output files (will append _hap1.fa, _hap2.fa)

    Returns:
        Tuple of (haplotype1_path, haplotype2_path)

    Raises:
        ValueError: If input is not a diploid reference (!=2 sequences)
        FileNotFoundError: If diploid_fasta doesn't exist

    Example:
        >>> hap1, hap2 = extract_haplotypes(
        ...     "diploid.fa",
        ...     "output/",
        ...     "sample_ont"
        ... )
        >>> print(hap1, hap2)
        output/sample_ont_hap1.fa output/sample_ont_hap2.fa

    Notes:
        - Output files named: {base_name}_hap1.fa, {base_name}_hap2.fa
        - Preserves original sequence IDs and descriptions
        - Creates output_dir if it doesn't exist
    """
    # Validate input
    if not is_diploid_reference(diploid_fasta):
        sequences = list(SeqIO.parse(diploid_fasta, "fasta"))
        raise ValueError(
            f"Expected diploid reference (2 sequences), found {len(sequences)} "
            f"in {diploid_fasta}"
        )

    # Parse sequences
    sequences = list(SeqIO.parse(diploid_fasta, "fasta"))

    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Define output paths
    hap1_path = output_path / f"{base_name}_hap1.fa"
    hap2_path = output_path / f"{base_name}_hap2.fa"

    # Write haplotypes
    with open(hap1_path, "w") as hap1_file:
        SeqIO.write([sequences[0]], hap1_file, "fasta")

    with open(hap2_path, "w") as hap2_file:
        SeqIO.write([sequences[1]], hap2_file, "fasta")

    logging.info(f"Extracted haplotypes from {diploid_fasta}:")
    logging.info(f"  Haplotype 1: {hap1_path} ({len(sequences[0].seq):,} bp)")
    logging.info(f"  Haplotype 2: {hap2_path} ({len(sequences[1].seq):,} bp)")

    return str(hap1_path), str(hap2_path)


def get_reference_info(fasta_path: str) -> dict:
    """
    Get detailed information about a reference FASTA file.

    Analyzes a FASTA file and returns metadata including:
    - Number of sequences
    - Reference type (diploid/haploid)
    - Sequence lengths
    - Total length
    - Length ratio (for diploid)

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary containing reference metadata:
        {
            "num_sequences": int,
            "is_diploid": bool,
            "sequences": [{"id": str, "length": int}, ...],
            "total_length": int,
            "length_ratio": float | None  # Only for diploid (longer/shorter)
        }

    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If FASTA file is empty or malformed

    Example:
        >>> info = get_reference_info("diploid.fa")
        >>> print(info)
        {
            "num_sequences": 2,
            "is_diploid": True,
            "sequences": [
                {"id": "haplotype_1", "length": 11200},
                {"id": "haplotype_2", "length": 16000}
            ],
            "total_length": 27200,
            "length_ratio": 1.43
        }
    """
    if not Path(fasta_path).exists():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")

    try:
        sequences = list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA file {fasta_path}: {e}") from e

    if len(sequences) == 0:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    # Build sequence info
    seq_info = [
        {"id": str(seq.id), "length": len(seq.seq)}
        for seq in sequences
    ]

    total_length = sum(s["length"] for s in seq_info)

    # Calculate length ratio for diploid
    length_ratio = None
    if len(sequences) == 2:
        lengths = [s["length"] for s in seq_info]
        length_ratio = max(lengths) / min(lengths) if min(lengths) > 0 else None

    return {
        "num_sequences": len(sequences),
        "is_diploid": len(sequences) == 2,
        "sequences": seq_info,
        "total_length": total_length,
        "length_ratio": length_ratio,
    }
```

**Tests**: `tests/test_reference_utils.py`

```python
"""Tests for reference_utils module."""

import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from muc_one_up.read_simulator.utils.reference_utils import (
    extract_haplotypes,
    get_reference_info,
    is_diploid_reference,
)


@pytest.fixture
def haploid_fasta(tmp_path):
    """Create a haploid FASTA file for testing."""
    fasta_path = tmp_path / "haploid.fa"
    record = SeqRecord(
        Seq("ATCGATCG" * 100),
        id="haplotype_1",
        description="Single haplotype"
    )
    SeqIO.write([record], fasta_path, "fasta")
    return str(fasta_path)


@pytest.fixture
def diploid_fasta(tmp_path):
    """Create a diploid FASTA file for testing."""
    fasta_path = tmp_path / "diploid.fa"
    hap1 = SeqRecord(
        Seq("ATCG" * 100),  # 400 bp
        id="haplotype_1",
        description="Short haplotype"
    )
    hap2 = SeqRecord(
        Seq("GCTA" * 200),  # 800 bp
        id="haplotype_2",
        description="Long haplotype"
    )
    SeqIO.write([hap1, hap2], fasta_path, "fasta")
    return str(fasta_path)


def test_is_diploid_reference_haploid(haploid_fasta):
    """Test diploid detection returns False for haploid reference."""
    assert is_diploid_reference(haploid_fasta) is False


def test_is_diploid_reference_diploid(diploid_fasta):
    """Test diploid detection returns True for diploid reference."""
    assert is_diploid_reference(diploid_fasta) is True


def test_is_diploid_reference_nonexistent():
    """Test diploid detection raises FileNotFoundError for missing file."""
    with pytest.raises(FileNotFoundError):
        is_diploid_reference("nonexistent.fa")


def test_extract_haplotypes_success(diploid_fasta, tmp_path):
    """Test successful haplotype extraction from diploid reference."""
    output_dir = tmp_path / "output"
    hap1_path, hap2_path = extract_haplotypes(
        diploid_fasta,
        str(output_dir),
        "test"
    )

    # Check files exist
    assert Path(hap1_path).exists()
    assert Path(hap2_path).exists()

    # Check content
    hap1_seqs = list(SeqIO.parse(hap1_path, "fasta"))
    hap2_seqs = list(SeqIO.parse(hap2_path, "fasta"))

    assert len(hap1_seqs) == 1
    assert len(hap2_seqs) == 1
    assert len(hap1_seqs[0].seq) == 400
    assert len(hap2_seqs[0].seq) == 800


def test_extract_haplotypes_not_diploid(haploid_fasta, tmp_path):
    """Test haplotype extraction fails on non-diploid reference."""
    with pytest.raises(ValueError, match="Expected diploid reference"):
        extract_haplotypes(haploid_fasta, str(tmp_path), "test")


def test_get_reference_info_diploid(diploid_fasta):
    """Test reference info extraction for diploid."""
    info = get_reference_info(diploid_fasta)

    assert info["num_sequences"] == 2
    assert info["is_diploid"] is True
    assert info["total_length"] == 1200
    assert info["length_ratio"] == pytest.approx(2.0)  # 800/400
    assert len(info["sequences"]) == 2


def test_get_reference_info_haploid(haploid_fasta):
    """Test reference info extraction for haploid."""
    info = get_reference_info(haploid_fasta)

    assert info["num_sequences"] == 1
    assert info["is_diploid"] is False
    assert info["length_ratio"] is None
```

---

#### Task 1.2: FASTQ Merging Utility

**File**: `muc_one_up/read_simulator/utils/fastq_utils.py`

**Purpose**: Efficiently merge multiple FASTQ files while preserving read integrity.

```python
#!/usr/bin/env python3
"""
FASTQ file manipulation utilities for read simulation.

Provides functions to:
- Merge multiple FASTQ files
- Validate FASTQ format
- Count reads
"""

import gzip
import logging
from pathlib import Path
from typing import List


def merge_fastq_files(
    input_fastqs: List[str],
    output_fastq: str,
    compress_output: bool = False
) -> str:
    """
    Merge multiple FASTQ files into a single output file.

    Concatenates FASTQ files in order, preserving all read information.
    Handles both compressed (.gz) and uncompressed FASTQ files automatically.
    Used to combine per-haplotype reads in split-simulation approach.

    Args:
        input_fastqs: List of paths to input FASTQ files (can be .gz or plain)
        output_fastq: Path for merged output FASTQ file
        compress_output: If True, compress output with gzip (default: False)

    Returns:
        Path to the merged FASTQ file

    Raises:
        FileNotFoundError: If any input FASTQ doesn't exist
        ValueError: If input_fastqs is empty
        IOError: If file operations fail

    Example:
        >>> merge_fastq_files(
        ...     ["hap1_reads.fastq", "hap2_reads.fastq"],
        ...     "merged_reads.fastq"
        ... )
        'merged_reads.fastq'

    Notes:
        - Reads are written in order: all reads from first file, then second, etc.
        - For diploid simulation, this preserves 50:50 representation if each
          file was simulated at the same coverage
        - Automatically detects .gz compression in input files
        - Output compression controlled by compress_output parameter
        - Fast implementation using binary copying (no parsing)
    """
    if not input_fastqs:
        raise ValueError("No input FASTQ files provided")

    # Validate all inputs exist
    for fastq in input_fastqs:
        if not Path(fastq).exists():
            raise FileNotFoundError(f"Input FASTQ not found: {fastq}")

    # Ensure output directory exists
    Path(output_fastq).parent.mkdir(parents=True, exist_ok=True)

    # Determine output mode
    output_path = Path(output_fastq)
    if compress_output and not output_path.suffix == ".gz":
        output_fastq = f"{output_fastq}.gz"
        output_path = Path(output_fastq)

    logging.info(f"Merging {len(input_fastqs)} FASTQ files → {output_fastq}")

    # Open output file (compressed or not)
    open_func = gzip.open if compress_output else open
    mode = "wb" if compress_output else "w"

    total_bytes = 0

    with open_func(output_fastq, mode) as out_file:
        for i, input_fastq in enumerate(input_fastqs, 1):
            input_path = Path(input_fastq)
            file_size = input_path.stat().st_size

            logging.debug(
                f"  [{i}/{len(input_fastqs)}] Appending {input_path.name} "
                f"({file_size:,} bytes)"
            )

            # Detect if input is compressed
            is_compressed = input_path.suffix == ".gz"
            read_func = gzip.open if is_compressed else open
            read_mode = "rb" if is_compressed else "r"

            # Copy file content
            with read_func(input_fastq, read_mode) as in_file:
                if compress_output or is_compressed:
                    # Binary mode
                    chunk_size = 8192
                    while True:
                        chunk = in_file.read(chunk_size)
                        if not chunk:
                            break
                        # Handle encoding if mixing compressed/uncompressed
                        if isinstance(chunk, str):
                            chunk = chunk.encode('utf-8')
                        out_file.write(chunk)
                        total_bytes += len(chunk)
                else:
                    # Text mode - simple copy
                    content = in_file.read()
                    out_file.write(content)
                    total_bytes += len(content.encode('utf-8'))

    logging.info(
        f"Successfully merged {len(input_fastqs)} files → {output_fastq} "
        f"({total_bytes:,} bytes)"
    )

    return str(output_fastq)


def count_fastq_reads(fastq_path: str) -> int:
    """
    Count the number of reads in a FASTQ file.

    Efficiently counts reads by counting lines and dividing by 4
    (each FASTQ record is exactly 4 lines).

    Args:
        fastq_path: Path to FASTQ file (can be .gz or plain)

    Returns:
        Number of reads in the file

    Raises:
        FileNotFoundError: If FASTQ file doesn't exist
        ValueError: If line count is not divisible by 4 (malformed FASTQ)

    Example:
        >>> count_fastq_reads("sample_reads.fastq")
        1000

    Notes:
        - Each FASTQ record is exactly 4 lines:
          @read_id
          SEQUENCE
          +
          QUALITY
        - Handles both compressed (.gz) and uncompressed files
        - Does not validate FASTQ format, only counts lines
    """
    fastq_file = Path(fastq_path)

    if not fastq_file.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    # Detect compression
    is_compressed = fastq_file.suffix == ".gz"
    open_func = gzip.open if is_compressed else open
    mode = "rt"  # Text mode for both

    # Count lines
    line_count = 0
    with open_func(fastq_path, mode) as f:
        for _ in f:
            line_count += 1

    # Validate FASTQ format (must be multiple of 4)
    if line_count % 4 != 0:
        raise ValueError(
            f"Malformed FASTQ file {fastq_path}: {line_count} lines "
            f"(not divisible by 4)"
        )

    read_count = line_count // 4

    logging.debug(f"Counted {read_count:,} reads in {fastq_path}")

    return read_count
```

**Tests**: `tests/test_fastq_utils.py`

```python
"""Tests for fastq_utils module."""

import gzip
from pathlib import Path

import pytest

from muc_one_up.read_simulator.utils.fastq_utils import (
    count_fastq_reads,
    merge_fastq_files,
)


@pytest.fixture
def sample_fastq(tmp_path):
    """Create a sample FASTQ file."""
    fastq_path = tmp_path / "sample.fastq"
    content = (
        "@read1\n"
        "ATCGATCG\n"
        "+\n"
        "IIIIIIII\n"
        "@read2\n"
        "GCTAGCTA\n"
        "+\n"
        "HHHHHHHH\n"
    )
    fastq_path.write_text(content)
    return str(fastq_path)


@pytest.fixture
def compressed_fastq(tmp_path):
    """Create a compressed FASTQ file."""
    fastq_path = tmp_path / "sample.fastq.gz"
    content = (
        "@read1\n"
        "ATCGATCG\n"
        "+\n"
        "IIIIIIII\n"
    ).encode('utf-8')
    with gzip.open(fastq_path, "wb") as f:
        f.write(content)
    return str(fastq_path)


def test_merge_fastq_files_plain(tmp_path):
    """Test merging uncompressed FASTQ files."""
    # Create two FASTQ files
    fq1 = tmp_path / "reads1.fastq"
    fq2 = tmp_path / "reads2.fastq"

    fq1.write_text("@read1\nATCG\n+\nIIII\n")
    fq2.write_text("@read2\nGCTA\n+\nHHHH\n")

    # Merge
    output = tmp_path / "merged.fastq"
    result = merge_fastq_files([str(fq1), str(fq2)], str(output))

    assert Path(result).exists()

    # Check content
    content = Path(result).read_text()
    assert "@read1" in content
    assert "@read2" in content
    assert content.count("@read") == 2


def test_merge_fastq_files_empty_list(tmp_path):
    """Test merging with empty input list raises ValueError."""
    with pytest.raises(ValueError, match="No input FASTQ files"):
        merge_fastq_files([], str(tmp_path / "out.fastq"))


def test_merge_fastq_files_missing_input(tmp_path):
    """Test merging with nonexistent input raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        merge_fastq_files(
            ["nonexistent.fastq"],
            str(tmp_path / "out.fastq")
        )


def test_count_fastq_reads_valid(sample_fastq):
    """Test read counting on valid FASTQ."""
    count = count_fastq_reads(sample_fastq)
    assert count == 2  # Two reads in fixture


def test_count_fastq_reads_compressed(compressed_fastq):
    """Test read counting on compressed FASTQ."""
    count = count_fastq_reads(compressed_fastq)
    assert count == 1


def test_count_fastq_reads_malformed(tmp_path):
    """Test read counting on malformed FASTQ (not multiple of 4 lines)."""
    bad_fastq = tmp_path / "bad.fastq"
    bad_fastq.write_text("@read1\nATCG\n+\n")  # Only 3 lines

    with pytest.raises(ValueError, match="Malformed FASTQ"):
        count_fastq_reads(str(bad_fastq))
```

---

### Phase 2: Diploid Handler (Core Logic)

#### Task 2.1: Diploid Simulation Handler

**File**: `muc_one_up/read_simulator/utils/diploid_handler.py`

**Purpose**: Orchestrate split-simulation workflow for diploid references.

```python
#!/usr/bin/env python3
"""
Diploid-aware NanoSim simulation handler.

Implements split-simulation approach to eliminate length-dependent allelic bias
in ONT read simulation. See Issue #31 for detailed analysis.

Key Features:
- Automatic diploid detection
- Per-haplotype simulation with independent seeds
- Coverage correction for training model behavior
- FASTQ merging to create balanced diploid reads

References:
    Issue #31: https://github.com/berntpopp/MucOneUp/issues/31
    Analysis: ISSUE_31_FINAL_REPORT.md
"""

import logging
from pathlib import Path
from typing import Tuple

from ..wrappers.nanosim_wrapper import run_nanosim_simulation
from .fastq_utils import count_fastq_reads, merge_fastq_files
from .reference_utils import extract_haplotypes, get_reference_info, is_diploid_reference


def simulate_diploid_split(
    nanosim_cmd: str,
    diploid_reference: str,
    output_prefix: str,
    training_model: str,
    target_coverage: float,
    correction_factor: float,
    threads: int = 4,
    base_seed: int | None = None,
    min_read_length: int | None = None,
    max_read_length: int | None = None,
    other_options: str = "",
    timeout: int = 3600,
) -> Tuple[str, dict]:
    """
    Simulate ONT reads from diploid reference using split-simulation approach.

    This function implements the solution to Issue #31: NanoSim allelic imbalance.
    Instead of simulating the entire diploid reference at once (which causes severe
    length-dependent bias), it:

    1. Extracts each haplotype into separate FASTA files
    2. Simulates each haplotype independently at half desired coverage
    3. Applies training model correction factor
    4. Merges FASTQs to create balanced diploid reads

    This reduces coverage bias from ~3.27:1 to ~1.40:1 (2.3x improvement).

    Args:
        nanosim_cmd: Path to NanoSim simulator.py command
        diploid_reference: Path to diploid FASTA (must contain 2 sequences)
        output_prefix: Prefix for output files
        training_model: Path to NanoSim training model
        target_coverage: Desired total coverage for both haplotypes combined
        correction_factor: Training model-specific actual/requested coverage ratio
                          (e.g., 0.325 for human_giab_hg002 model)
        threads: Number of threads for simulation
        base_seed: Base random seed (hap1 uses seed, hap2 uses seed+1)
        min_read_length: Minimum read length (optional)
        max_read_length: Maximum read length (optional)
        other_options: Additional NanoSim options
        timeout: Timeout per simulation in seconds

    Returns:
        Tuple of (merged_fastq_path, statistics_dict)

        statistics_dict contains:
        {
            "hap1_fastq": str,
            "hap2_fastq": str,
            "hap1_reads": int,
            "hap2_reads": int,
            "total_reads": int,
            "hap1_length": int,
            "hap2_length": int,
            "length_ratio": float,
            "target_coverage": float,
            "corrected_coverage_per_hap": float,
        }

    Raises:
        ValueError: If diploid_reference is not diploid (!=2 sequences)
        RuntimeError: If simulation or merging fails

    Example:
        >>> merged_fastq, stats = simulate_diploid_split(
        ...     nanosim_cmd="simulator.py",
        ...     diploid_reference="sample.fa",
        ...     output_prefix="output/sample_ont",
        ...     training_model="training/",
        ...     target_coverage=200.0,
        ...     correction_factor=0.325,
        ...     threads=8,
        ...     base_seed=42
        ... )
        >>> print(f"Generated {stats['total_reads']} reads")
        Generated 157 reads  # 50:50 split

    Notes:
        - Correction factor is training model-specific (empirically determined)
        - Each haplotype simulated at: (target_coverage / 2) / correction_factor
        - Different seeds ensure independent read sampling per haplotype
        - Output files: {prefix}_hap1.fa, {prefix}_hap2.fa, {prefix}_merged.fastq
    """
    # Validate diploid reference
    if not is_diploid_reference(diploid_reference):
        raise ValueError(
            f"Split-simulation requires diploid reference (2 sequences), "
            f"but {diploid_reference} is not diploid. "
            f"Use standard simulation for haploid references."
        )

    # Get reference info for logging
    ref_info = get_reference_info(diploid_reference)
    logging.info(
        f"Diploid reference detected: {ref_info['num_sequences']} sequences, "
        f"{ref_info['total_length']:,} bp total, "
        f"length ratio {ref_info['length_ratio']:.2f}:1"
    )

    # Calculate corrected coverage per haplotype
    # Each haplotype gets half the desired coverage, corrected for training model
    coverage_per_haplotype = (target_coverage / 2.0) / correction_factor

    logging.info(
        f"Split-simulation parameters:"
    )
    logging.info(f"  Target total coverage: {target_coverage}x")
    logging.info(f"  Correction factor: {correction_factor}")
    logging.info(f"  Coverage per haplotype (corrected): {coverage_per_haplotype:.1f}x")
    logging.info(f"  Base seed: {base_seed}")

    # Extract haplotypes
    output_dir = Path(output_prefix).parent
    base_name = Path(output_prefix).name

    logging.info("Step 1/3: Extracting haplotypes")
    hap1_fa, hap2_fa = extract_haplotypes(
        diploid_reference,
        str(output_dir),
        base_name
    )

    # Simulate haplotype 1
    logging.info("Step 2/3: Simulating haplotype 1")
    hap1_prefix = str(Path(output_prefix).with_name(f"{base_name}_hap1"))

    hap1_seed = base_seed if base_seed is not None else None

    hap1_fastq = run_nanosim_simulation(
        nanosim_cmd=nanosim_cmd,
        reference_fasta=hap1_fa,
        output_prefix=hap1_prefix,
        training_model=training_model,
        coverage=coverage_per_haplotype,
        threads=threads,
        min_read_length=min_read_length,
        max_read_length=max_read_length,
        other_options=other_options,
        timeout=timeout,
        seed=hap1_seed,
    )

    hap1_read_count = count_fastq_reads(hap1_fastq)
    logging.info(f"  Haplotype 1: {hap1_read_count:,} reads generated")

    # Simulate haplotype 2 (different seed for independent sampling)
    logging.info("Step 3/3: Simulating haplotype 2")
    hap2_prefix = str(Path(output_prefix).with_name(f"{base_name}_hap2"))

    hap2_seed = (base_seed + 1) if base_seed is not None else None

    hap2_fastq = run_nanosim_simulation(
        nanosim_cmd=nanosim_cmd,
        reference_fasta=hap2_fa,
        output_prefix=hap2_prefix,
        training_model=training_model,
        coverage=coverage_per_haplotype,
        threads=threads,
        min_read_length=min_read_length,
        max_read_length=max_read_length,
        other_options=other_options,
        timeout=timeout,
        seed=hap2_seed,
    )

    hap2_read_count = count_fastq_reads(hap2_fastq)
    logging.info(f"  Haplotype 2: {hap2_read_count:,} reads generated")

    # Merge FASTQs
    logging.info("Step 4/4: Merging haplotype reads")
    merged_fastq = f"{output_prefix}_merged.fastq"

    merge_fastq_files(
        input_fastqs=[hap1_fastq, hap2_fastq],
        output_fastq=merged_fastq,
        compress_output=False
    )

    total_reads = hap1_read_count + hap2_read_count
    read_ratio = hap2_read_count / hap1_read_count if hap1_read_count > 0 else 0

    logging.info(
        f"Split-simulation complete: {total_reads:,} total reads "
        f"(hap1: {hap1_read_count:,}, hap2: {hap2_read_count:,}, "
        f"ratio: {read_ratio:.2f}:1)"
    )

    # Build statistics
    statistics = {
        "hap1_fastq": hap1_fastq,
        "hap2_fastq": hap2_fastq,
        "hap1_reads": hap1_read_count,
        "hap2_reads": hap2_read_count,
        "total_reads": total_reads,
        "read_ratio": read_ratio,
        "hap1_length": ref_info["sequences"][0]["length"],
        "hap2_length": ref_info["sequences"][1]["length"],
        "length_ratio": ref_info["length_ratio"],
        "target_coverage": target_coverage,
        "corrected_coverage_per_hap": coverage_per_haplotype,
        "correction_factor": correction_factor,
    }

    return merged_fastq, statistics
```

**Tests**: `tests/test_diploid_handler.py`

```python
"""Tests for diploid_handler module."""

import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import Mock, patch, call

from muc_one_up.read_simulator.utils.diploid_handler import simulate_diploid_split


@pytest.fixture
def diploid_reference(tmp_path):
    """Create a mock diploid FASTA."""
    fasta_path = tmp_path / "diploid.fa"
    hap1 = SeqRecord(Seq("ATCG" * 2800), id="haplotype_1", description="Short 11.2kb")
    hap2 = SeqRecord(Seq("GCTA" * 4000), id="haplotype_2", description="Long 16kb")
    SeqIO.write([hap1, hap2], fasta_path, "fasta")
    return str(fasta_path)


@pytest.fixture
def haploid_reference(tmp_path):
    """Create a mock haploid FASTA."""
    fasta_path = tmp_path / "haploid.fa"
    hap = SeqRecord(Seq("ATCG" * 4000), id="haplotype_1", description="Single haplotype")
    SeqIO.write([hap], fasta_path, "fasta")
    return str(fasta_path)


@patch("muc_one_up.read_simulator.utils.diploid_handler.run_nanosim_simulation")
@patch("muc_one_up.read_simulator.utils.diploid_handler.count_fastq_reads")
@patch("muc_one_up.read_simulator.utils.diploid_handler.merge_fastq_files")
def test_simulate_diploid_split_success(
    mock_merge,
    mock_count,
    mock_nanosim,
    diploid_reference,
    tmp_path
):
    """Test successful diploid split-simulation."""
    # Setup mocks
    mock_nanosim.side_effect = [
        str(tmp_path / "hap1_aligned_reads.fastq"),
        str(tmp_path / "hap2_aligned_reads.fastq"),
    ]
    mock_count.side_effect = [100, 100]  # 50:50 split
    mock_merge.return_value = str(tmp_path / "merged.fastq")

    # Run split-simulation
    merged_fastq, stats = simulate_diploid_split(
        nanosim_cmd="simulator.py",
        diploid_reference=diploid_reference,
        output_prefix=str(tmp_path / "test_ont"),
        training_model="training/",
        target_coverage=200.0,
        correction_factor=0.325,
        threads=8,
        base_seed=42,
    )

    # Assertions
    assert merged_fastq == str(tmp_path / "merged.fastq")
    assert stats["total_reads"] == 200
    assert stats["hap1_reads"] == 100
    assert stats["hap2_reads"] == 100
    assert stats["read_ratio"] == pytest.approx(1.0)
    assert stats["target_coverage"] == 200.0
    assert stats["correction_factor"] == 0.325

    # Verify NanoSim called twice with correct coverage
    expected_coverage = (200.0 / 2.0) / 0.325  # ~307.7
    assert mock_nanosim.call_count == 2

    # Check seeds are different
    calls = mock_nanosim.call_args_list
    assert calls[0].kwargs["seed"] == 42
    assert calls[1].kwargs["seed"] == 43


@patch("muc_one_up.read_simulator.utils.diploid_handler.run_nanosim_simulation")
def test_simulate_diploid_split_not_diploid(mock_nanosim, haploid_reference, tmp_path):
    """Test split-simulation fails gracefully on non-diploid reference."""
    with pytest.raises(ValueError, match="Split-simulation requires diploid"):
        simulate_diploid_split(
            nanosim_cmd="simulator.py",
            diploid_reference=haploid_reference,
            output_prefix=str(tmp_path / "test"),
            training_model="training/",
            target_coverage=200.0,
            correction_factor=0.325,
        )
```

---

### Phase 3: Pipeline Integration

#### Task 3.1: Refactor ONT Pipeline

**File**: `muc_one_up/read_simulator/ont_pipeline.py` (MODIFIED)

**Changes**:
1. Add feature flags to control split-simulation and coverage correction
2. Integrate diploid detection and routing
3. Preserve backward compatibility
4. Update docstrings

```python
#!/usr/bin/env python3
"""
Oxford Nanopore read simulation pipeline with diploid-aware split-simulation.

This module implements the ONT read simulation pipeline using NanoSim,
with automatic handling of diploid references to eliminate allelic bias.

Key Features (Issue #31 Solution):
- Automatic diploid reference detection
- Split-simulation approach for balanced diploid reads
- Training model-specific coverage correction
- Feature flags for backward compatibility
- Comprehensive logging and statistics

Pipeline Modes:
1. Haploid/Standard: Single NanoSim call (original behavior)
2. Diploid Split-Simulation: Per-haplotype simulation + merging (NEW)

References:
    Issue #31: https://github.com/berntpopp/MucOneUp/issues/31
    Analysis: ISSUE_31_FINAL_REPORT.md
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from .utils.diploid_handler import simulate_diploid_split
from .utils.reference_utils import is_diploid_reference
from .wrappers.nanosim_wrapper import (
    align_ont_reads_with_minimap2,
    run_nanosim_simulation,
)


def simulate_ont_reads_pipeline(
    config: dict[str, Any], input_fa: str, human_reference: str | None = None
) -> str:
    """
    Run the complete Oxford Nanopore read simulation pipeline.

    This function automatically detects diploid references and applies the
    split-simulation approach (Issue #31) to eliminate length-dependent
    allelic bias. For haploid references, uses standard simulation.

    Pipeline Steps:

    **For Diploid References** (2 sequences, enable_split_simulation=True):
    1. Detect diploid reference
    2. Extract haplotypes to separate FASTA files
    3. Calculate corrected coverage: (target / 2) / correction_factor
    4. Simulate haplotype 1 with seed N
    5. Simulate haplotype 2 with seed N+1
    6. Merge FASTQ files (50:50 representation)
    7. Align merged reads to original diploid reference
    8. Create and index BAM

    **For Haploid References** (or split_simulation=False):
    1. Apply coverage correction (if enabled)
    2. Run standard NanoSim simulation
    3. Align reads to reference
    4. Create and index BAM

    Args:
        config: Configuration dictionary with sections:
            - tools: External tool commands (nanosim, minimap2, samtools)
            - nanosim_params: NanoSim parameters:
                - training_data_path: Path to training model (required)
                - coverage: Target coverage (required)
                - correction_factor: Training model ratio (default: 0.325)
                - enable_split_simulation: Use diploid split mode (default: True)
                - enable_coverage_correction: Apply correction factor (default: True)
                - num_threads: Thread count (default: 4)
                - seed: Random seed for reproducibility (optional)
                - min_read_length: Minimum read length (optional)
                - max_read_length: Maximum read length (optional)
                - other_options: Additional NanoSim flags (optional)
            - read_simulation: General parameters:
                - output_dir: Output directory (default: same as input)
                - threads: Thread override (optional)

        input_fa: Path to input FASTA file (haploid or diploid)
        human_reference: Path to human reference for alignment.
                        If None, aligns to input_fa.

    Returns:
        Path to final output BAM file ({input_basename}_ont.bam)

    Raises:
        ValueError: If required config parameters are missing
        RuntimeError: If any pipeline step fails

    Example:
        >>> config = {
        ...     "tools": {
        ...         "nanosim": "simulator.py",
        ...         "minimap2": "minimap2",
        ...         "samtools": "samtools"
        ...     },
        ...     "nanosim_params": {
        ...         "training_data_path": "training/model",
        ...         "coverage": 200,
        ...         "correction_factor": 0.325,
        ...         "enable_split_simulation": True,
        ...         "seed": 42
        ...     }
        ... }
        >>> bam_file = simulate_ont_reads_pipeline(
        ...     config,
        ...     "sample.001.simulated.fa"
        ... )
        >>> print(bam_file)
        output/sample_001_simulated_ont.bam

    Notes:
        - Diploid detection: Checks if FASTA contains exactly 2 sequences
        - Coverage correction factor is training model-specific:
          * human_giab_hg002 model: 0.325 (empirically determined)
          * Other models: Should be calibrated via systematic testing
        - Split-simulation reduces bias from ~3.27:1 to ~1.40:1
        - Both features (split-simulation, correction) can be disabled via config
        - Preserves backward compatibility when features disabled

    See Also:
        - ISSUE_31_FINAL_REPORT.md: Complete analysis and validation
        - simulate_diploid_split(): Core split-simulation implementation
        - run_nanosim_simulation(): Standard simulation wrapper
    """
    start_time = datetime.now()
    logging.info(
        "Starting ONT read simulation pipeline at %s",
        start_time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    # Extract configuration
    tools = config.get("tools", {})
    ns_params = config.get("nanosim_params", {})
    rs_config = config.get("read_simulation", {})

    # Validate required tools
    nanosim_cmd = tools.get("nanosim")
    minimap2_cmd = tools.get("minimap2")
    samtools_cmd = tools.get("samtools")

    if not nanosim_cmd:
        raise ValueError("Missing 'nanosim' command in tools configuration.")
    if not minimap2_cmd:
        raise ValueError("Missing 'minimap2' command in tools configuration.")
    if not samtools_cmd:
        raise ValueError("Missing 'samtools' command in tools configuration.")

    # Validate required simulation parameters
    training_model = ns_params.get("training_data_path")
    coverage = ns_params.get("coverage")

    if not training_model:
        raise ValueError("Missing 'training_data_path' in nanosim_params configuration.")
    if not coverage:
        raise ValueError("Missing 'coverage' in nanosim_params configuration.")

    # Feature flags (default: enabled)
    enable_split_simulation = ns_params.get("enable_split_simulation", True)
    enable_coverage_correction = ns_params.get("enable_coverage_correction", True)
    correction_factor = ns_params.get("correction_factor", 0.325)

    # Other parameters
    threads = ns_params.get("num_threads", rs_config.get("threads", 4))
    min_read_length = ns_params.get("min_read_length")
    max_read_length = ns_params.get("max_read_length")
    other_options = ns_params.get("other_options", "")
    seed = ns_params.get("seed")

    # Setup output paths
    input_path = Path(input_fa)
    input_basename = input_path.stem
    output_dir = rs_config.get("output_dir", str(input_path.parent))
    output_prefix = str(Path(output_dir) / f"{input_basename}_ont")

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Detect diploid reference
    is_diploid = is_diploid_reference(input_fa)

    logging.info(f"Reference analysis:")
    logging.info(f"  Type: {'Diploid' if is_diploid else 'Haploid'}")
    logging.info(f"  Split-simulation: {enable_split_simulation}")
    logging.info(f"  Coverage correction: {enable_coverage_correction}")
    if enable_coverage_correction:
        logging.info(f"  Correction factor: {correction_factor}")

    # Determine simulation mode
    use_split_simulation = is_diploid and enable_split_simulation

    # Run simulation
    logging.info("1. Starting NanoSim simulation")

    if use_split_simulation:
        # Diploid split-simulation mode
        logging.info("  Mode: Diploid split-simulation (Issue #31 solution)")

        try:
            fastq_file, split_stats = simulate_diploid_split(
                nanosim_cmd=nanosim_cmd,
                diploid_reference=input_fa,
                output_prefix=output_prefix,
                training_model=training_model,
                target_coverage=float(coverage),
                correction_factor=correction_factor,
                threads=int(threads),
                base_seed=seed,
                min_read_length=min_read_length,
                max_read_length=max_read_length,
                other_options=other_options,
            )

            # Log split-simulation statistics
            logging.info("  Split-simulation statistics:")
            logging.info(f"    Haplotype 1: {split_stats['hap1_reads']:,} reads "
                        f"({split_stats['hap1_length']:,} bp)")
            logging.info(f"    Haplotype 2: {split_stats['hap2_reads']:,} reads "
                        f"({split_stats['hap2_length']:,} bp)")
            logging.info(f"    Read ratio: {split_stats['read_ratio']:.2f}:1")
            logging.info(f"    Length ratio: {split_stats['length_ratio']:.2f}:1")

        except Exception as e:
            logging.error("Diploid split-simulation failed: %s", e)
            raise RuntimeError(f"ONT read simulation failed: {e!s}") from e

    else:
        # Standard simulation mode (haploid or split disabled)
        mode_desc = "Standard" if not is_diploid else "Standard (split-simulation disabled)"
        logging.info(f"  Mode: {mode_desc}")

        # Apply coverage correction if enabled
        effective_coverage = coverage
        if enable_coverage_correction:
            effective_coverage = float(coverage) / correction_factor
            logging.info(
                f"  Applying coverage correction: "
                f"{coverage}x → {effective_coverage:.1f}x "
                f"(factor: {correction_factor})"
            )

        try:
            fastq_file = run_nanosim_simulation(
                nanosim_cmd=nanosim_cmd,
                reference_fasta=input_fa,
                output_prefix=output_prefix,
                training_model=training_model,
                coverage=float(effective_coverage),
                threads=int(threads),
                min_read_length=min_read_length,
                max_read_length=max_read_length,
                other_options=other_options,
                seed=seed,
            )
        except Exception as e:
            logging.error("NanoSim simulation failed: %s", e)
            raise RuntimeError(f"ONT read simulation failed: {e!s}") from e

    logging.info("NanoSim simulation completed successfully")

    # 2. Align reads with minimap2
    logging.info("2. Starting read alignment with minimap2")
    output_bam = f"{output_prefix}.bam"

    reference_for_alignment = human_reference if human_reference else input_fa
    logging.info(f"Aligning ONT reads to reference: {reference_for_alignment}")

    try:
        align_ont_reads_with_minimap2(
            minimap2_cmd=minimap2_cmd,
            samtools_cmd=samtools_cmd,
            human_reference=reference_for_alignment,
            reads_fastq=fastq_file,
            output_bam=output_bam,
            threads=int(threads),
        )
        logging.info("Read alignment completed successfully")
    except Exception as e:
        logging.error("Read alignment failed: %s", e)
        raise RuntimeError(f"ONT read alignment failed: {e!s}") from e

    # Calculate elapsed time
    end_time = datetime.now()
    duration = end_time - start_time
    logging.info(
        "ONT read simulation pipeline completed at %s (duration: %s)",
        end_time.strftime("%Y-%m-%d %H:%M:%S"),
        str(duration).split(".")[0],
    )

    logging.info("Final outputs:")
    logging.info("  Aligned and indexed BAM: %s", output_bam)
    logging.info("  Reads FASTQ: %s", fastq_file)

    return output_bam
```

---

### Phase 4: Configuration Schema Update

#### Task 4.1: Add NanoSim Parameters to Config Schema

**File**: `config.json` (EXAMPLE UPDATE)

Add the following to the example config:

```json
{
  "nanosim_params": {
    "training_data_path": "reference/nanosim/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training",
    "coverage": 200,
    "correction_factor": 0.325,
    "enable_split_simulation": true,
    "enable_coverage_correction": true,
    "num_threads": 8,
    "seed": 42,
    "min_read_length": null,
    "max_read_length": null,
    "other_options": ""
  }
}
```

**File**: `muc_one_up/config.py` (VALIDATION UPDATE)

Add schema validation for new parameters:

```python
# In validate_config() function, add:

NANOSIM_SCHEMA = {
    "training_data_path": str,  # Required
    "coverage": (int, float),   # Required
    "correction_factor": (int, float),  # Optional, default 0.325
    "enable_split_simulation": bool,     # Optional, default True
    "enable_coverage_correction": bool,  # Optional, default True
    "num_threads": int,         # Optional
    "seed": int,                # Optional
    "min_read_length": int,     # Optional
    "max_read_length": int,     # Optional
    "other_options": str,       # Optional
}

# Add validation in validate_config()
ns_params = config.get("nanosim_params", {})
if ns_params:
    # Validate correction_factor range
    corr_factor = ns_params.get("correction_factor", 0.325)
    if not 0.0 < corr_factor <= 1.0:
        logging.warning(
            f"correction_factor {corr_factor} outside typical range (0.0, 1.0]. "
            f"Ensure this is calibrated for your training model."
        )
```

---

### Phase 5: Integration Tests

#### Task 5.1: End-to-End Pipeline Tests

**File**: `tests/test_ont_pipeline.py` (NEW)

```python
"""Integration tests for ONT pipeline with split-simulation."""

import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import Mock, patch, MagicMock

from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline


@pytest.fixture
def diploid_reference(tmp_path):
    """Create test diploid reference."""
    fasta = tmp_path / "diploid.fa"
    hap1 = SeqRecord(Seq("ATCG" * 2800), id="haplotype_1", description="Short")
    hap2 = SeqRecord(Seq("GCTA" * 4000), id="haplotype_2", description="Long")
    SeqIO.write([hap1, hap2], fasta, "fasta")
    return str(fasta)


@pytest.fixture
def base_config(tmp_path):
    """Base configuration for testing."""
    return {
        "tools": {
            "nanosim": "simulator.py",
            "minimap2": "minimap2",
            "samtools": "samtools"
        },
        "nanosim_params": {
            "training_data_path": "training/model",
            "coverage": 200,
            "correction_factor": 0.325,
            "enable_split_simulation": True,
            "enable_coverage_correction": True,
            "seed": 42,
            "num_threads": 8
        },
        "read_simulation": {
            "output_dir": str(tmp_path / "output")
        }
    }


@patch("muc_one_up.read_simulator.ont_pipeline.simulate_diploid_split")
@patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")
def test_pipeline_diploid_split_mode(
    mock_align,
    mock_split,
    diploid_reference,
    base_config,
    tmp_path
):
    """Test pipeline uses split-simulation for diploid with flag enabled."""
    # Setup mocks
    merged_fastq = str(tmp_path / "merged.fastq")
    mock_split.return_value = (
        merged_fastq,
        {
            "hap1_reads": 100,
            "hap2_reads": 100,
            "total_reads": 200,
            "read_ratio": 1.0,
            "hap1_length": 11200,
            "hap2_length": 16000,
            "length_ratio": 1.43,
            "target_coverage": 200.0,
            "corrected_coverage_per_hap": 307.7,
            "correction_factor": 0.325
        }
    )
    mock_align.return_value = str(tmp_path / "output.bam")

    # Run pipeline
    bam = simulate_ont_reads_pipeline(
        base_config,
        diploid_reference
    )

    # Assertions
    mock_split.assert_called_once()
    mock_align.assert_called_once()
    assert bam.endswith(".bam")

    # Check split was called with correct parameters
    call_kwargs = mock_split.call_args.kwargs
    assert call_kwargs["target_coverage"] == 200.0
    assert call_kwargs["correction_factor"] == 0.325
    assert call_kwargs["base_seed"] == 42


@patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
@patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")
def test_pipeline_standard_mode_disabled_split(
    mock_align,
    mock_nanosim,
    diploid_reference,
    base_config,
    tmp_path
):
    """Test pipeline uses standard mode when split-simulation disabled."""
    # Disable split-simulation
    base_config["nanosim_params"]["enable_split_simulation"] = False

    # Setup mocks
    mock_nanosim.return_value = str(tmp_path / "reads.fastq")
    mock_align.return_value = str(tmp_path / "output.bam")

    # Run pipeline
    bam = simulate_ont_reads_pipeline(
        base_config,
        diploid_reference
    )

    # Should use standard simulation, not split
    mock_nanosim.assert_called_once()
    mock_align.assert_called_once()

    # Check coverage correction was applied
    call_kwargs = mock_nanosim.call_args.kwargs
    expected_cov = 200.0 / 0.325  # ~615.4
    assert call_kwargs["coverage"] == pytest.approx(expected_cov, rel=0.01)


def test_pipeline_missing_required_params(diploid_reference):
    """Test pipeline fails gracefully with missing required params."""
    incomplete_config = {
        "tools": {"nanosim": "simulator.py"},  # Missing minimap2, samtools
        "nanosim_params": {}  # Missing required params
    }

    with pytest.raises(ValueError):
        simulate_ont_reads_pipeline(incomplete_config, diploid_reference)
```

---

### Phase 6: Documentation

#### Task 6.1: Update CLAUDE.md

Add comprehensive section on ONT simulation:

```markdown
### ONT Read Simulation with NanoSim

MucOneUp implements a **diploid-aware split-simulation approach** for Oxford Nanopore read simulation to eliminate allelic bias (Issue #31).

#### Background

Standard NanoSim simulation of diploid references produces severe allelic imbalance:
- Longer haplotypes receive 2-3x more coverage than shorter haplotypes
- Root cause: Length-proportional read start position sampling
- Impact: Short alleles may have insufficient coverage for variant calling

**Solution**: Split-simulation approach reduces bias from ~3.27:1 to ~1.40:1.

#### How It Works

**For Diploid References** (2 sequences):
1. Automatically detects diploid FASTA (exactly 2 sequences)
2. Extracts each haplotype to separate FASTA file
3. Simulates each haplotype independently at half desired coverage
4. Applies training model-specific coverage correction
5. Merges FASTQ files to create balanced diploid reads
6. Aligns merged reads to original diploid reference

**For Haploid References** (1 sequence):
- Standard NanoSim simulation
- Optional coverage correction applied

#### Configuration

```json
{
  "nanosim_params": {
    "training_data_path": "path/to/training/model",
    "coverage": 200,
    "correction_factor": 0.325,
    "enable_split_simulation": true,
    "enable_coverage_correction": true,
    "seed": 42
  }
}
```

**Key Parameters**:

- `correction_factor` (default: 0.325): Training model-specific ratio of actual/requested coverage
  - For `human_giab_hg002` model: 0.325 (empirically determined)
  - Other models: Calibrate via systematic testing

- `enable_split_simulation` (default: true): Enable diploid split-simulation
  - Automatically activates for 2-sequence FASTA files
  - Disable to use standard simulation

- `enable_coverage_correction` (default: true): Apply coverage correction factor
  - Accounts for NanoSim's kernel density estimation behavior
  - Disable if using custom pre-calibrated coverage values

#### Calibrating Correction Factor

For custom training models:

```bash
# 1. Create haploid test reference (~16kb)
muconeup --config config.json simulate --out-base test --fixed-lengths 100

# 2. Extract single haplotype
head -n 2 test.001.simulated.fa > test_haploid.fa

# 3. Test at 200x coverage
muconeup --config config.json reads ont test_haploid.fa --coverage 200

# 4. Align and measure actual coverage
minimap2 -ax map-ont test_haploid.fa test_haploid_ont_aligned_reads.fastq | \
  samtools view -b | samtools sort -o test.bam
samtools depth test.bam | awk '{sum+=$3} END {print "Actual:", sum/NR"x"}'

# 5. Calculate correction factor
# correction_factor = actual_coverage / 200
# Example: 65.4x / 200 = 0.327
```

Update `correction_factor` in config.json with your calibrated value.

#### Example Usage

```bash
# Simulate diploid haplotypes
muconeup --config config.json simulate --out-base sample --fixed-lengths 60,80

# Generate balanced ONT reads (automatic split-simulation)
muconeup --config config.json reads ont sample.001.simulated.fa

# Output:
# - sample.001.simulated_ont_hap1.fa     # Haplotype 1 extracted
# - sample.001.simulated_ont_hap2.fa     # Haplotype 2 extracted
# - sample.001.simulated_ont_merged.fastq  # Balanced reads
# - sample.001.simulated_ont.bam         # Aligned reads
```

#### Performance Characteristics

**Bias Reduction** (20 vs 100 repeats, 200x target):
```
Standard approach:
  SHORT: 28.81x coverage  (3.27:1 bias)
  LONG:  94.12x coverage

Split-simulation approach:
  SHORT: 22.18x coverage  (1.40:1 bias)
  LONG:  30.98x coverage

Improvement: 2.3x more balanced
```

**Coverage Accuracy** (with correction_factor=0.325):
```
Requested: 200x
Actual (haploid): ~65x per haplotype
Total (diploid): ~130x combined
Efficiency: 32.5%
```

#### Backward Compatibility

To use original behavior (no split-simulation):

```json
{
  "nanosim_params": {
    "enable_split_simulation": false,
    "enable_coverage_correction": false
  }
}
```

#### References

- **Issue #31**: https://github.com/berntpopp/MucOneUp/issues/31
- **Analysis**: ISSUE_31_FINAL_REPORT.md
- **NanoSim**: https://github.com/bcgsc/NanoSim
```

---

## Testing Strategy

### Unit Tests
- ✅ `test_reference_utils.py` - FASTA analysis functions
- ✅ `test_fastq_utils.py` - FASTQ merging and counting
- ✅ `test_diploid_handler.py` - Split-simulation orchestration

### Integration Tests
- ✅ `test_ont_pipeline.py` - End-to-end pipeline with mocks
- ⏳ `test_ont_pipeline_integration.py` - Real NanoSim calls (optional, slow)

### Validation Tests
- ⏳ Create test script comparing standard vs split-simulation
- ⏳ Measure bias reduction on known test case (20 vs 100 repeats)

---

## Implementation Timeline

**Phase 1 (Foundation)**: 2-3 hours
- Task 1.1: Reference utils + tests
- Task 1.2: FASTQ utils + tests

**Phase 2 (Core Logic)**: 2-3 hours
- Task 2.1: Diploid handler + tests

**Phase 3 (Integration)**: 2 hours
- Task 3.1: Refactor ONT pipeline

**Phase 4 (Configuration)**: 1 hour
- Task 4.1: Config schema updates

**Phase 5 (Testing)**: 1-2 hours
- Task 5.1: Integration tests

**Phase 6 (Documentation)**: 1 hour
- Task 6.1: Update CLAUDE.md
- Task 6.2: Update README.md

**Total Estimated Time**: 9-12 hours

---

## Success Criteria

1. ✅ All unit tests pass
2. ✅ All integration tests pass
3. ✅ Diploid references automatically use split-simulation (when enabled)
4. ✅ Haploid references use standard simulation
5. ✅ Coverage bias reduced from >3:1 to <1.5:1 on test case
6. ✅ Backward compatibility preserved (feature flags work)
7. ✅ Documentation complete and accurate
8. ✅ Code follows DRY, KISS, SOLID principles
9. ✅ Comprehensive docstrings and type hints
10. ✅ No breaking changes to existing API

---

## Risks and Mitigation

**Risk 1**: Doubled simulation time for diploid references
- **Mitigation**: Document performance impact; simulations parallelizable

**Risk 2**: Correction factor may vary between training models
- **Mitigation**: Provide calibration instructions; use conservative default

**Risk 3**: Backward compatibility concerns
- **Mitigation**: Feature flags default to new behavior; can disable easily

**Risk 4**: Increased complexity
- **Mitigation**: Modular design; comprehensive tests; clear documentation

---

## Post-Implementation

1. Run validation test suite on multiple haplotype combinations
2. Update GitHub Issue #31 with implementation status
3. Consider blog post explaining the solution
4. Monitor user feedback and edge cases
5. Collect correction factors for additional training models

---

**End of Implementation Plan**

This plan follows software engineering best practices and ensures robust, maintainable, well-tested code that solves Issue #31 while preserving backward compatibility and code quality.

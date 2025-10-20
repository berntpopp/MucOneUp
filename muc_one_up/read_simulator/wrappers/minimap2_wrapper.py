#!/usr/bin/env python3
"""
Wrapper module for minimap2 alignment operations.

This module provides a generic minimap2 wrapper with parameterized preset
support, enabling alignment of ONT, PacBio CLR, and PacBio HiFi reads.

Key Functions:
    align_reads_with_minimap2: Generic minimap2 alignment with configurable preset

Key Features:
    - Parameterized preset support (map-ont, map-pb, map-hifi)
    - Automatic SAM→BAM conversion via samtools
    - BAM sorting and indexing
    - Thread control for performance
    - Timeout handling for large datasets

Example:
    Align PacBio HiFi reads::

        from muc_one_up.read_simulator.wrappers.minimap2_wrapper import align_reads_with_minimap2
        from muc_one_up.read_simulator.constants import MINIMAP2_PRESET_PACBIO_HIFI

        bam_file = align_reads_with_minimap2(
            minimap2_cmd="minimap2",
            samtools_cmd="samtools",
            reference="hg38.fa",
            reads_fastq="hifi_reads.fastq",
            output_bam="aligned.bam",
            preset=MINIMAP2_PRESET_PACBIO_HIFI,
            threads=8
        )

Notes:
    - Replaces hardcoded ONT-specific alignment function
    - Eliminates code duplication across pipelines (DRY principle)
    - Preset parameter enables technology-specific optimization
    - All intermediate SAM files are automatically cleaned up

References:
    - minimap2 presets: https://github.com/lh3/minimap2#usage
    - SAM/BAM format: https://samtools.github.io/hts-specs/SAMv1.pdf
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..constants import (
    DEFAULT_ALIGNMENT_TIMEOUT,
    DEFAULT_MINIMAP2_THREADS,
    MINIMAP2_PRESET_ONT,
)
from ..utils import run_command
from .samtools_wrapper import convert_sam_to_bam, sort_and_index_bam


def align_reads_with_minimap2(
    minimap2_cmd: str,
    samtools_cmd: str,
    reference: str,
    reads_fastq: str,
    output_bam: str,
    preset: str = MINIMAP2_PRESET_ONT,
    threads: int = DEFAULT_MINIMAP2_THREADS,
    timeout: int = DEFAULT_ALIGNMENT_TIMEOUT,
) -> str:
    """
    Align long reads to reference genome using minimap2 with configurable preset.

    This generic wrapper supports multiple sequencing technologies (ONT, PacBio CLR,
    PacBio HiFi) via parameterized preset selection. Automatically converts SAM to
    sorted, indexed BAM format.

    Args:
        minimap2_cmd: Path to the minimap2 executable.
        samtools_cmd: Path to the samtools executable.
        reference: Reference genome FASTA file path.
        reads_fastq: Input FASTQ file containing reads to align.
        output_bam: Output BAM file path for aligned reads.
        preset: Minimap2 alignment preset (default: map-ont).
                Valid presets:
                - map-ont: Oxford Nanopore reads
                - map-pb: PacBio CLR reads
                - map-hifi: PacBio HiFi reads (CCS)
        threads: Number of threads to use (default: 4).
        timeout: Timeout in seconds (default: 3600).

    Returns:
        Path to the sorted, indexed output BAM file.

    Raises:
        ExternalToolError: If minimap2 or samtools command fails.
        FileOperationError: If input files don't exist or output creation fails.

    Example:
        Align ONT reads (default preset)::

            from muc_one_up.read_simulator.wrappers.minimap2_wrapper import align_reads_with_minimap2

            bam_file = align_reads_with_minimap2(
                minimap2_cmd="minimap2",
                samtools_cmd="samtools",
                reference="hg38.fa",
                reads_fastq="ont_reads.fastq",
                output_bam="aligned.bam"
            )

        Align PacBio HiFi reads::

            from muc_one_up.read_simulator.constants import MINIMAP2_PRESET_PACBIO_HIFI

            bam_file = align_reads_with_minimap2(
                minimap2_cmd="minimap2",
                samtools_cmd="samtools",
                reference="hg38.fa",
                reads_fastq="hifi_reads.fastq",
                output_bam="aligned.bam",
                preset=MINIMAP2_PRESET_PACBIO_HIFI,
                threads=16
            )

    Notes:
        - Input FASTQ can be gzip-compressed (.fastq.gz)
        - Reference FASTA is automatically indexed if .mmi doesn't exist
        - Intermediate SAM file is cleaned up after conversion
        - Output BAM is coordinate-sorted and indexed
        - Preset parameter enables technology-specific optimization:
          * map-ont: Higher indel error model for ONT
          * map-pb: Optimized for PacBio CLR error profile
          * map-hifi: Low error model for high-accuracy HiFi reads
    """
    # Validate input files exist
    reference_path = Path(reference)
    if not reference_path.exists():
        raise FileOperationError(f"Reference genome not found: {reference}")

    reads_fastq_path = Path(reads_fastq)
    if not reads_fastq_path.exists():
        raise FileOperationError(f"Input FASTQ file not found: {reads_fastq}")

    # Create intermediate SAM filename
    output_sam = str(Path(output_bam).with_suffix(".sam"))

    # Run minimap2 alignment with specified preset
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        minimap2_cmd,
        "-ax",
        preset,  # Technology-specific preset (map-ont, map-pb, map-hifi)
        "-t",
        threads,  # Threads (build_tool_command handles conversion)
        reference,
        reads_fastq,
        "-o",
        output_sam,
    )

    logging.info(f"Aligning reads with minimap2 (preset: {preset}): {reads_fastq} → {output_sam}")
    logging.info(f"Using reference: {reference}")
    run_command(cmd, timeout=timeout, stderr_prefix="[minimap2] ", stderr_log_level=logging.INFO)

    # Validate SAM output exists and is non-empty
    output_sam_path = Path(output_sam)
    if not output_sam_path.exists() or output_sam_path.stat().st_size == 0:
        raise FileOperationError(
            f"minimap2 alignment failed: Output SAM {output_sam} missing or empty"
        )

    # Convert SAM to BAM using reusable samtools wrapper
    logging.info(f"Converting SAM to BAM: {output_sam} → {output_bam}")
    convert_sam_to_bam(
        samtools_cmd=samtools_cmd,
        input_sam=output_sam,
        output_bam=output_bam,
        threads=threads,
        timeout=timeout,
    )

    # Sort and index the BAM file using reusable samtools wrapper
    logging.info(f"Sorting and indexing BAM: {output_bam}")
    sort_and_index_bam(
        samtools_exe=samtools_cmd,
        input_bam=output_bam,
        output_bam=output_bam,
        threads=threads,
    )

    # Clean up intermediate SAM file
    try:
        if output_sam_path.exists():
            output_sam_path.unlink()
            logging.debug(f"Removed intermediate SAM file: {output_sam}")
    except Exception as e:
        logging.warning(f"Could not remove intermediate SAM file {output_sam}: {e}")

    # Validate final BAM exists and is non-empty
    output_bam_path = Path(output_bam)
    if not output_bam_path.exists() or output_bam_path.stat().st_size == 0:
        raise FileOperationError(
            f"Alignment pipeline failed: Output BAM {output_bam} missing or empty"
        )

    # Validate BAM index exists
    bam_index = Path(f"{output_bam}.bai")
    if not bam_index.exists():
        logging.warning(f"BAM index not found: {bam_index}")

    logging.info(f"Alignment complete: {output_bam}")
    return output_bam


def validate_preset(preset: str) -> None:
    """
    Validate minimap2 preset string.

    Args:
        preset: Minimap2 alignment preset to validate.

    Raises:
        ValueError: If preset is not one of the valid options.

    Example:
        Validate before alignment::

            from muc_one_up.read_simulator.wrappers.minimap2_wrapper import validate_preset
            from muc_one_up.read_simulator.constants import MINIMAP2_PRESET_PACBIO_HIFI

            validate_preset(MINIMAP2_PRESET_PACBIO_HIFI)  # OK
            validate_preset("invalid-preset")  # Raises ValueError

    Notes:
        - Valid presets imported from constants module
        - Used for early validation in pipeline functions
        - Prevents silent failures from typos in preset names
    """
    from ..constants import (
        MINIMAP2_PRESET_ONT,
        MINIMAP2_PRESET_PACBIO_CLR,
        MINIMAP2_PRESET_PACBIO_HIFI,
    )

    valid_presets = {
        MINIMAP2_PRESET_ONT,
        MINIMAP2_PRESET_PACBIO_CLR,
        MINIMAP2_PRESET_PACBIO_HIFI,
    }

    if preset not in valid_presets:
        valid_list = ", ".join(sorted(valid_presets))
        raise ValueError(f"Invalid minimap2 preset: '{preset}'. Valid options: {valid_list}")

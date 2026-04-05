"""Illumina fragment preparation stage (pipeline stages 1-8).

Extracts the first 8 stages of the Illumina read simulation pipeline into a
standalone function with a clean, testable interface.

Stages covered:
    1. Replace Ns in input FASTA (reseq replaceN)
    2. Validate reseq_model config
    3. Generate systematic errors (reseq illuminaPE)
    4. Convert FASTA to 2bit (faToTwoBit)
    5. Resolve sample BAM from assembly context or config
    6. Extract subset reference from sample BAM (samtools)
    7. Align subset reference with pblat
    8. Simulate fragments (Wessim2-style ported logic)
    9. Create reads from fragments (reseq seqToIllumina)
   10. Split interleaved FASTQ into R1/R2 pairs
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ...exceptions import ConfigurationError
from ..assembly_context import AssemblyContext
from ..fragment_simulation import simulate_fragments
from ..wrappers.reseq_wrapper import (
    create_reads,
    generate_systematic_errors,
    replace_Ns,
    split_reads,
)
from ..wrappers.samtools_wrapper import extract_subset_reference
from ..wrappers.ucsc_tools_wrapper import fa_to_twobit, run_pblat
from . import FragmentResult

logger = logging.getLogger(__name__)


def prepare_fragments(
    tools: dict[str, str],
    rs_config: dict[str, Any],
    input_fa: str,
    output_dir: Path,
    output_base: str,
    assembly_ctx: AssemblyContext,
    source_tracker: Any | None = None,
) -> FragmentResult:
    """Run Illumina fragment preparation (pipeline stages 1-8).

    Produces paired-end FASTQ files ready for downstream alignment.

    Args:
        tools: Mapping of tool names to executable paths/commands.
        rs_config: Read-simulation configuration dictionary.
        input_fa: Path to the input FASTA file.
        output_dir: Directory where output and intermediate files are written.
        output_base: Base name used for all output and intermediate file names.
        assembly_ctx: Resolved assembly context (provides sample_bam etc.).
        source_tracker: Optional read source tracker; when not None, fragment
            origins are written to ``{output_base}_fragment_origins.tsv``.

    Returns:
        FragmentResult with R1/R2 FASTQ paths and a list of intermediate files.

    Raises:
        ConfigurationError: If ``reseq_model`` is missing from *rs_config*, or
            no sample BAM can be resolved from *assembly_ctx* or *rs_config*.
    """
    intermediate_files: list[str] = []

    # Stage 1: Replace Ns
    no_ns_fa = str(Path(output_dir) / f"_{output_base}_noNs.fa")
    logger.info("1. Replacing Ns in FASTA")
    replace_Ns(input_fa, no_ns_fa, tools)
    intermediate_files.append(no_ns_fa)

    # Stage 2: Validate reseq_model
    reseq_model = rs_config.get("reseq_model")
    if not reseq_model:
        raise ConfigurationError("reseq_model not specified in config 'read_simulation' section")

    # Stage 3: Generate systematic errors
    syser_fq = str(Path(output_dir) / f"_{output_base}_syser.fq")
    logger.info("2. Generating systematic errors")
    generate_systematic_errors(no_ns_fa, reseq_model, syser_fq, tools)
    intermediate_files.append(syser_fq)

    # Stage 4: Convert FASTA to 2bit
    twobit_file = str(Path(output_dir) / f"_{output_base}_noNs.2bit")
    logger.info("3. Converting FASTA to 2bit format")
    fa_to_twobit(no_ns_fa, twobit_file, tools)
    intermediate_files.append(twobit_file)

    # Stage 5: Resolve sample BAM
    sample_bam = assembly_ctx.sample_bam
    if not sample_bam:
        sample_bam = rs_config.get("sample_bam")
        if sample_bam:
            logger.warning(
                "Assembly-specific sample BAM not found for %s. Using generic sample_bam: %s",
                assembly_ctx.assembly_name,
                sample_bam,
            )

    if not sample_bam:
        raise ConfigurationError(
            f"No sample BAM specified for {assembly_ctx.assembly_name}. "
            f"Add 'sample_bam_{assembly_ctx.assembly_name}' or 'sample_bam' to config"
        )

    # Stage 6: Extract subset reference from BAM
    subset_ref = str(Path(output_dir) / f"_{output_base}_subset_ref.fa")
    logger.info("4. Extracting subset reference from BAM")
    collated_bam = extract_subset_reference(sample_bam, subset_ref, tools)
    intermediate_files.append(subset_ref)
    intermediate_files.append(collated_bam)

    # Stage 7: Run pblat alignment
    psl_file = str(Path(output_dir) / f"_{output_base}_alignment.psl")
    threads = rs_config.get("threads", 4)
    pblat_threads = min(rs_config.get("pblat_threads", 24), threads)
    logger.info("5. Running pblat alignment")
    run_pblat(
        twobit_file,
        subset_ref,
        psl_file,
        tools,
        threads=pblat_threads,
        min_score=rs_config.get("pblat_min_score", 95),
        min_identity=rs_config.get("pblat_min_identity", 95),
    )

    # Stage 8: Simulate fragments (Wessim2-style ported logic)
    fragments_fa = str(Path(output_dir) / f"_{output_base}_fragments.fa")
    read_number = rs_config.get("read_number", 100000)
    fragment_size = rs_config.get("fragment_size", 350)
    fragment_sd = rs_config.get("fragment_sd", 50)
    min_fragment = rs_config.get("min_fragment", 200)
    bind = rs_config.get("binding_min", 0.5)
    seed = rs_config.get("seed")
    logger.info("6. Simulating fragments (Wessim2-style, ported logic)")
    fragment_origins_path = (
        str(Path(output_dir) / f"{output_base}_fragment_origins.tsv")
        if source_tracker is not None
        else None
    )
    simulate_fragments(
        no_ns_fa,
        syser_fq,
        psl_file,
        read_number,
        fragment_size,
        fragment_sd,
        min_fragment,
        bind,
        fragments_fa,
        seed=seed,
        fragment_origins_path=fragment_origins_path,
    )
    intermediate_files.append(psl_file)
    intermediate_files.append(fragments_fa)

    # Stage 9: Create reads from fragments
    reads_fq = str(Path(output_dir) / f"_{output_base}_reads.fq")
    logger.info("7. Creating reads from fragments")
    create_reads(
        fragments_fa,
        reseq_model,
        reads_fq,
        threads,
        tools,
    )

    # Stage 10: Split reads into paired FASTQ files
    reads_fq1 = str(Path(output_dir) / f"{output_base}_R1.fastq.gz")
    reads_fq2 = str(Path(output_dir) / f"{output_base}_R2.fastq.gz")
    logger.info("8. Splitting reads into paired FASTQ files")
    split_reads(reads_fq, reads_fq1, reads_fq2)
    intermediate_files.append(reads_fq)

    return FragmentResult(
        r1_fastq=reads_fq1,
        r2_fastq=reads_fq2,
        intermediate_files=intermediate_files,
    )

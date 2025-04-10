#!/usr/bin/env python3
"""
Read simulation pipeline orchestrator.

This module is the core implementation of the MUC1 read simulation pipeline.
It orchestrates the multi-step process of simulating Illumina reads from a
FASTA file by coordinating calls to wrapper modules and fragment simulation.

Key features:
- Timeout handling for all external tool calls to prevent pipeline hanging
- Consistent file naming with clear distinction between outputs and intermediates
- Comprehensive error handling and validation
- Optional downsampling to achieve target coverage

Implementation details:
- All intermediate files are prefixed with underscore (_) to prevent accidental deletion
- All external tools are run with appropriate timeouts (60-300s depending on complexity)
- The final outputs use the input filename as base for consistency

For usage information, see the main read_simulation.py module.
"""

import logging
import os
import sys
from datetime import datetime
from typing import Dict, Any

# Import wrapper modules
from .wrappers.reseq_wrapper import (
    replace_Ns,
    generate_systematic_errors,
    create_reads,
    split_reads,
)
from .wrappers.samtools_wrapper import (
    extract_subset_reference,
    calculate_vntr_coverage,
    calculate_target_coverage,
    downsample_bam,
    downsample_entire_bam,
)
from .wrappers.bwa_wrapper import align_reads
from .wrappers.ucsc_tools_wrapper import fa_to_twobit, run_pblat

# Import fragment simulation
from .fragment_simulation import simulate_fragments

# Import utilities
from .utils import check_external_tools, cleanup_files


def simulate_reads_pipeline(config: Dict[str, Any], input_fa: str) -> str:
    """
    Run the complete read simulation pipeline.

    Pipeline steps:
    1. Replace Ns in the simulated FASTA using reseq replaceN
    2. Generate systematic errors using reseq illuminaPE
    3. Convert FASTA to 2bit format using faToTwoBit
    4. Extract a subset reference from a sample BAM using samtools
    5. Align the 2bit file to the subset reference using pblat
    6. Simulate fragments using a ported version of w‑Wessim2
    7. Create reads from fragments using reseq seqToIllumina
    8. Split the interleaved FASTQ into paired FASTQ files
    9. Align the reads to a human reference using BWA MEM
    10. Optionally downsample to a target coverage
    11. Clean up intermediate files

    Args:
        config: Dictionary containing "tools" and "read_simulation" sections.
              Required config parameters:
              - tools: Dictionary of tool commands (reseq, faToTwoBit, samtools, pblat, bwa)
              - read_simulation: Dictionary of simulation parameters
                - reseq_model: Path to reseq model file
                - sample_bam: Sample BAM file to extract reference
                - human_reference: Reference genome FASTA
              Optional parameters:
              - seqtoillumina_timeout: Custom timeout for seqToIllumina tool (default: 120s)
              - output_dir: Output directory (default: input file directory)
              - read_number: Number of read pairs to simulate (default: 100000)
              - fragment_size: Mean fragment size (default: 350)
              - fragment_sd: Fragment size standard deviation (default: 50)
              - min_fragment: Minimum fragment size (default: 200)
        input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).

    Returns:
        Path to the final output BAM file ({input_basename}.bam).

    Raises:
        SystemExit: If any step in the pipeline fails.
    """
    # Start time
    start_time = datetime.now()
    logging.info(
        "Starting read simulation pipeline at %s",
        start_time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    # Extract tool commands and read simulation config
    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})

    # Validate external tools
    check_external_tools(tools)

    # Create output directory if it doesn't exist
    output_dir = rs_config.get("output_dir", os.path.dirname(input_fa))
    os.makedirs(output_dir, exist_ok=True)

    # Define the common base name for all output files
    output_base = os.path.basename(input_fa).replace(".fa", "").replace(".fasta", "")

    # Use consistent naming for all output files based on the output_base
    reads_fq1 = rs_config.get(
        "output_fastq1", os.path.join(output_dir, f"{output_base}_R1.fastq.gz")
    )
    reads_fq2 = rs_config.get(
        "output_fastq2", os.path.join(output_dir, f"{output_base}_R2.fastq.gz")
    )
    output_bam = rs_config.get(
        "output_bam", os.path.join(output_dir, f"{output_base}.bam")
    )

    # Log the output filenames for clarity
    logging.info("Output filenames:")
    logging.info("  BAM: %s", output_bam)
    logging.info("  FASTQ pair: %s, %s", reads_fq1, reads_fq2)

    # Stage 1: Replace Ns in the input FASTA
    # Define all intermediate filenames consistently with underscore prefix
    # to distinguish them from final outputs
    noNs_fa = os.path.join(output_dir, f"_{output_base}_noNs.fa")
    logging.info("1. Replacing Ns in FASTA")
    replace_Ns(input_fa, noNs_fa, tools)

    # Stage 2: Generate systematic errors
    reseq_model = rs_config.get("reseq_model")
    if not reseq_model:
        logging.error("Error: reseq_model not specified in config.")
        sys.exit(1)
    syser_fq = os.path.join(output_dir, f"_{output_base}_syser.fq")
    logging.info("2. Generating systematic errors")
    generate_systematic_errors(noNs_fa, reseq_model, syser_fq, tools)

    # Stage 3: Convert FASTA to 2bit
    twobit_file = os.path.join(output_dir, f"_{output_base}_noNs.2bit")
    logging.info("3. Converting FASTA to 2bit format")
    fa_to_twobit(noNs_fa, twobit_file, tools)

    # Stage 4: Extract subset reference from BAM
    # Get the reference assembly setting from the config
    reference_assembly = rs_config.get("reference_assembly", "hg38")

    # Get the appropriate sample_bam based on the reference assembly
    sample_bam_key = f"sample_bam_{reference_assembly}"
    sample_bam = rs_config.get(sample_bam_key)

    # Fall back to the generic sample_bam if assembly-specific one not found
    if not sample_bam:
        sample_bam = rs_config.get("sample_bam")
        if sample_bam:
            logging.warning(
                f"Assembly-specific sample BAM not found for {reference_assembly}. "
                f"Using generic sample_bam: {sample_bam}"
            )

    if not sample_bam:
        logging.error(
            f"Error: No sample BAM specified for {reference_assembly}. "
            f"Please add 'sample_bam_{reference_assembly}' or 'sample_bam' to the config."
        )
        sys.exit(1)

    logging.info(f"Using sample BAM for {reference_assembly}: {sample_bam}")
    subset_ref = os.path.join(output_dir, f"_{output_base}_subset_ref.fa")
    logging.info("4. Extracting subset reference from BAM")
    collated_bam = extract_subset_reference(sample_bam, subset_ref, tools)

    # Stage 5: Run pblat alignment
    psl_file = os.path.join(output_dir, f"_{output_base}_alignment.psl")
    threads = rs_config.get("threads", 4)
    pblat_threads = min(rs_config.get("pblat_threads", 24), threads)
    logging.info("5. Running pblat alignment")
    run_pblat(
        twobit_file,
        subset_ref,
        psl_file,
        tools,
        threads=pblat_threads,
        minScore=rs_config.get("pblat_min_score", 95),
        minIdentity=rs_config.get("pblat_min_identity", 95),
    )

    # Stage 6: Simulate fragments
    fragments_fa = os.path.join(output_dir, f"_{output_base}_fragments.fa")
    read_number = rs_config.get("read_number", 100000)
    fragment_size = rs_config.get("fragment_size", 350)
    fragment_sd = rs_config.get("fragment_sd", 50)
    min_fragment = rs_config.get("min_fragment", 200)
    bind = rs_config.get("binding_min", 0.5)
    logging.info("6. Simulating fragments (w-Wessim2)")
    simulate_fragments(
        noNs_fa,
        syser_fq,
        psl_file,
        read_number,
        fragment_size,
        fragment_sd,
        min_fragment,
        bind,
        fragments_fa,
    )

    # Stage 7: Create reads from fragments
    reads_fq = os.path.join(output_dir, f"_{output_base}_reads.fq")
    logging.info("7. Creating reads from fragments")
    # Use a timeout of 120 seconds for seqToIllumina - the tool tends to keep running
    # indefinitely even after producing complete output
    seqtoillumina_timeout = rs_config.get(
        "seqtoillumina_timeout", 120
    )  # Default 120 seconds
    create_reads(
        fragments_fa,
        reseq_model,
        reads_fq,
        threads,
        tools,
        timeout=seqtoillumina_timeout,
    )

    # Stage 8: Split reads into paired FASTQ files
    logging.info("8. Splitting reads into paired FASTQ files")
    split_reads(reads_fq, reads_fq1, reads_fq2)

    # Stage 9: Align reads
    human_reference = rs_config.get("human_reference")
    if not human_reference:
        logging.error("Error: human_reference not specified in config.")
        sys.exit(1)
    logging.info("9. Aligning reads to reference")
    align_reads(reads_fq1, reads_fq2, human_reference, output_bam, tools, threads)

    # Stage 10: Optionally downsample
    downsample_target = rs_config.get("downsample_target")
    if downsample_target:
        logging.info("10. Downsampling to target coverage")
        mode = rs_config.get("downsample_mode", "vntr").strip().lower()
        if mode == "vntr":
            reference_assembly = rs_config.get("reference_assembly", "hg38")
            vntr_region_key = f"vntr_region_{reference_assembly}"
            vntr_region = rs_config.get(vntr_region_key)
            if not vntr_region:
                logging.error(
                    f"VNTR region not specified in config for {reference_assembly}."
                )
                sys.exit(1)
            current_cov = calculate_vntr_coverage(
                tools["samtools"],
                output_bam,
                vntr_region,
                threads,
                os.path.dirname(output_bam),
                os.path.basename(output_bam).replace(".bam", ""),
            )
            region_info = vntr_region
        elif mode == "non_vntr":
            bed_file = rs_config.get("sample_target_bed")
            if not bed_file:
                logging.error(
                    "For non-VNTR downsampling, 'sample_target_bed' must be provided in config."
                )
                sys.exit(1)
            current_cov = calculate_target_coverage(
                tools["samtools"],
                output_bam,
                bed_file,
                threads,
                os.path.dirname(output_bam),
                os.path.basename(output_bam).replace(".bam", ""),
            )
            region_info = f"BED file: {bed_file}"
        else:
            logging.error(
                "Invalid downsample_mode in config; use 'vntr' or 'non_vntr'."
            )
            sys.exit(1)
        if current_cov > downsample_target:
            fraction = downsample_target / current_cov
            fraction = min(max(fraction, 0.0), 1.0)
            logging.info(
                "Downsampling BAM from %.2fx to target %dx (fraction: %.4f) based on %s",
                current_cov,
                downsample_target,
                fraction,
                region_info,
            )
            downsampled_bam = os.path.join(
                os.path.dirname(output_bam),
                os.path.basename(output_bam).replace(".bam", "_downsampled.bam"),
            )
            if mode == "vntr":
                downsample_bam(
                    tools["samtools"],
                    output_bam,
                    downsampled_bam,
                    vntr_region,
                    fraction,
                    rs_config.get("downsample_seed", 42),
                    threads,
                )
            else:
                downsample_entire_bam(
                    tools["samtools"],
                    output_bam,
                    downsampled_bam,
                    fraction,
                    rs_config.get("downsample_seed", 42),
                    threads,
                )
            output_bam = downsampled_bam
        else:
            logging.info(
                "Current coverage (%.2fx) is below the target; no downsampling performed.",
                current_cov,
            )
    else:
        logging.info("10. Skipping downsampling (no target coverage specified)")

    # Clean up intermediate files
    end_time = datetime.now()
    duration = end_time - start_time
    logging.info(
        "Read simulation pipeline completed at %s (duration: %s)",
        end_time.strftime("%Y-%m-%d %H:%M:%S"),
        str(duration).split(".")[0],  # Remove microseconds for cleaner output
    )
    logging.info("Final outputs:")
    logging.info("  Aligned and indexed BAM: %s", output_bam)
    logging.info("  Paired FASTQ files (gzipped): %s and %s", reads_fq1, reads_fq2)

    # Clean up intermediate files
    intermediates = [
        noNs_fa,
        syser_fq,
        twobit_file,
        subset_ref,
        psl_file,
        fragments_fa,
        reads_fq,
        collated_bam,
    ]

    # Double-check to make sure output files are not in the intermediate list
    if output_bam in intermediates:
        logging.warning("Prevented deletion of output BAM file: %s", output_bam)
        intermediates.remove(output_bam)

    cleanup_files(intermediates)

    return output_bam


# For direct command-line use (backward compatibility)
if __name__ == "__main__":
    import json

    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    if len(sys.argv) != 3:
        print("Usage: python pipeline.py <config.json> <input_fasta>")
        sys.exit(1)

    config_file = sys.argv[1]
    input_fa = sys.argv[2]

    with open(config_file, "r") as fh:
        config = json.load(fh)

    simulate_reads_pipeline(config, input_fa)

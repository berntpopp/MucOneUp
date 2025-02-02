#!/usr/bin/env python3
"""
bam_anonymizer.py

This script subsets a BAM file to a specific region, then reheaders and
anonymizes the resulting BAM using GATK's ReadAnonymizer tool. Finally, it
performs a reheader step to remove @pg and @rg lines from the output BAM.
The final anonymized and reheadered BAM is written with a name based on the
provided target design.

Example usage:
    python bam_anonymizer.py \
        --input-bam /path/to/input.bam \
        --target-design twist_v2 \
        --ref /path/to/reference.fasta \
        --keep-intermediates

Required external tools:
  - samtools (for subsetting, reheadering, and indexing BAM files)
  - gatk (for the ReadAnonymizer tool)

The default subsetting region is set to the hg38 MUC1 VNTR region:
    chr1:155184000-155194000

You may override the region via the --region argument.
"""

import argparse
import logging
import os
import subprocess
import sys


def subset_bam(input_bam: str, region: str, output_bam: str) -> None:
    """
    Subset the input BAM to a given region using samtools view and index it.

    Args:
        input_bam: Path to the input BAM file.
        region: Region to subset (e.g. "chr1:155184000-155194000").
        output_bam: Path for the output subset BAM.
    """
    logging.info("Subsetting BAM: %s to region: %s", input_bam, region)
    view_cmd = [
        "samtools",
        "view",
        "-b",  # output BAM format
        "-P",  # keep paired reads together
        input_bam,
        region,
        "-o",
        output_bam,
    ]
    logging.debug("Running command: %s", " ".join(view_cmd))
    subprocess.run(view_cmd, check=True)

    # Index the subset BAM file.
    index_cmd = ["samtools", "index", output_bam]
    logging.debug("Indexing subset BAM with command: %s", " ".join(index_cmd))
    subprocess.run(index_cmd, check=True)


def run_gatk_read_anonymizer(input_bam: str, output_bam: str, reference: str) -> None:
    """
    Run GATK's ReadAnonymizer to remove sensitive genetic information.

    Args:
        input_bam: Path to the input BAM file (typically the subset BAM).
        output_bam: Path for the output anonymized BAM file.
        reference: Path to the reference FASTA file.
    """
    logging.info("Running GATK ReadAnonymizer on: %s", input_bam)
    gatk_cmd = [
        "gatk",
        "ReadAnonymizer",
        "-I",
        input_bam,
        "-O",
        output_bam,
        "-R",
        reference,
    ]
    logging.debug("Running command: %s", " ".join(gatk_cmd))
    subprocess.run(gatk_cmd, check=True)


def run_samtools_reheader(input_bam: str, output_bam: str) -> None:
    """
    Reheader the BAM to remove @pg and @rg lines, then index it.

    Args:
        input_bam: Path to the input BAM file (e.g. output from GATK ReadAnonymizer).
        output_bam: Path for the reheadered output BAM file.
    """
    logging.info("Reheadering BAM: %s", input_bam)
    cmd_reheader = [
        "samtools", "reheader",
        "-P",
        "-c", "grep -v ^@PG | grep -v ^@RG | grep -v ^@CO",
        input_bam
    ]
    logging.debug("Reheader command: %s", " ".join(cmd_reheader))
    with open(output_bam, 'w') as fout:
        subprocess.run(cmd_reheader, stdout=fout, check=True)

    # Index the reheadered BAM file.
    index_cmd = ["samtools", "index", output_bam]
    logging.debug("Indexing reheadered BAM with command: %s", " ".join(index_cmd))
    subprocess.run(index_cmd, check=True)


def index_bam(bam_file: str) -> None:
    """
    Index the provided BAM file using samtools.

    Args:
        bam_file: Path to the BAM file to be indexed.
    """
    logging.info("Indexing final BAM: %s", bam_file)
    index_cmd = ["samtools", "index", bam_file]
    subprocess.run(index_cmd, check=True)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Subset a BAM file to a region and anonymize it using GATK ReadAnonymizer. "
            "Then reheader the anonymized BAM to remove @pg and @rg lines. "
            "The final output file is named based on the provided target design."
        )
    )
    parser.add_argument(
        "--input-bam",
        required=True,
        help="Path to the input BAM file.",
    )
    parser.add_argument(
        "--target-design",
        required=True,
        help=(
            "Target design name to be used for naming the output file "
            "(e.g. twist_v2)."
        ),
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Path to the reference FASTA file (required by GATK ReadAnonymizer).",
    )
    parser.add_argument(
        "--region",
        default="chr1:155184000-155194000",
        help=(
            "Genomic region to subset from the BAM file. "
            "Default is 'chr1:155184000-155194000' (hg38 MUC1 VNTR region)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to write output files. Default is current directory.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help=(
            "Keep intermediate files (e.g. the subset BAM and its index) for debugging purposes."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Set logging level. Default is INFO.",
    )
    return parser.parse_args()


def main() -> None:
    """Main function for the BAM anonymization pipeline."""
    args = parse_args()
    logging.basicConfig(
        level=args.log_level.upper(),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # Verify that the input BAM and reference exist.
    if not os.path.isfile(args.input_bam):
        logging.error("Input BAM file not found: %s", args.input_bam)
        sys.exit(1)
    if not os.path.isfile(args.ref):
        logging.error("Reference FASTA file not found: %s", args.ref)
        sys.exit(1)

    # Ensure the output directory exists.
    os.makedirs(args.output_dir, exist_ok=True)

    # Define file paths.
    subset_bam_path = os.path.join(args.output_dir, f"{args.target_design}_subset.bam")
    # This file is produced by GATK ReadAnonymizer (but will be reheadered next).
    final_bam_path = os.path.join(args.output_dir, f"{args.target_design}.bam")
    # Final output after reheadering.
    reheadered_bam_path = os.path.join(args.output_dir, f"{args.target_design}_reheader.bam")

    try:
        # 1. Subset the BAM file to the specified region.
        subset_bam(args.input_bam, args.region, subset_bam_path)

        # 2. Run GATK ReadAnonymizer on the subset BAM.
        run_gatk_read_anonymizer(subset_bam_path, final_bam_path, args.ref)

        # 3. Perform final reheadering on the anonymized BAM.
        run_samtools_reheader(final_bam_path, reheadered_bam_path)
        # Optionally, remove the non-reheadered file.
        os.remove(final_bam_path)
        final_bam_path = reheadered_bam_path

        logging.info("Anonymized and reheadered BAM file created successfully: %s", final_bam_path)
    except subprocess.CalledProcessError as err:
        logging.error("A subprocess failed: %s", err)
        sys.exit(1)
    finally:
        if args.keep_intermediates:
            logging.info("Keeping intermediate files (--keep-intermediates set): %s", subset_bam_path)
        else:
            # Remove the temporary subset BAM file and its index.
            if os.path.exists(subset_bam_path):
                try:
                    os.remove(subset_bam_path)
                    subset_index = subset_bam_path + ".bai"
                    if os.path.exists(subset_index):
                        os.remove(subset_index)
                    logging.debug("Removed intermediate files.")
                except Exception as err:
                    logging.warning("Could not remove temporary file(s): %s", err)


if __name__ == "__main__":
    main()

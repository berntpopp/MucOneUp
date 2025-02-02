#!/usr/bin/env python3
"""
bam_anonymizer.py

This script subsets a BAM file to a specific region, replaces the read group
information with generic values using GATK AddOrReplaceReadGroups, anonymizes
the resulting BAM using GATK's ReadAnonymizer, and reheaders the anonymized
BAM to remove unwanted header lines (such as @PG, @RG, and @CO). The final
output BAM is named exactly according to the provided target design (e.g.
twist_v2.bam).

After processing, a pseudonymisation output CSV is written with columns:
  old_name, old_md5, new_name, new_md5

Example usage:
    python bam_anonymizer.py \
        --input-bam /path/to/input.bam \
        --target-design twist_v2 \
        --ref /path/to/reference.fasta \
        --keep-intermediates

Required external tools:
  - samtools (for subsetting, reheadering, and indexing BAM files)
  - gatk (for the ReadAnonymizer and AddOrReplaceReadGroups tools)

The default subsetting region is set to the hg38 MUC1 VNTR region:
    chr1:155184000-155194000

You may override the region via the --region argument.
"""

import argparse
import csv
import hashlib
import logging
import os
import subprocess
import sys


def compute_md5(filepath: str, chunk_size: int = 1048576) -> str:
    """
    Compute the MD5 hex digest of the given file in a memory-efficient way.

    Args:
        filepath: Path to the file.
        chunk_size: Number of bytes to read at a time (default: 1MB).

    Returns:
        The MD5 hexadecimal digest as a string.
    """
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()


def subset_bam(input_bam: str, region: str, output_bam: str) -> None:
    """
    Subset the input BAM to a given region using samtools view and index it.

    Args:
        input_bam: Path to the input BAM file.
        region: Genomic region to subset (e.g. "chr1:155184000-155194000").
        output_bam: Path for the output subset BAM.
    """
    logging.info("Subsetting BAM: %s to region: %s", input_bam, region)
    view_cmd = [
        "samtools", "view", "-b", "-P", input_bam, region, "-o", output_bam
    ]
    logging.debug("Running command: %s", " ".join(view_cmd))
    subprocess.run(view_cmd, check=True)
    index_cmd = ["samtools", "index", output_bam]
    logging.debug("Indexing subset BAM with command: %s", " ".join(index_cmd))
    subprocess.run(index_cmd, check=True)


def run_gatk_add_or_replace_read_groups(input_bam: str, output_bam: str) -> None:
    """
    Replace the read groups in the BAM with generic values using GATK AddOrReplaceReadGroups.

    Args:
        input_bam: Path to the input BAM (e.g. the subset BAM).
        output_bam: Path for the output BAM with generic read group information.
    """
    logging.info("Replacing read groups on: %s", input_bam)
    cmd = [
        "gatk", "AddOrReplaceReadGroups",
        "-I", input_bam,
        "-O", output_bam,
        "-RGID", "4",
        "-RGLB", "lib1",
        "-RGPL", "ILLUMINA",
        "-RGPU", "unit1",
        "-RGSM", "20"
    ]
    logging.debug("Running command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_gatk_read_anonymizer(input_bam: str, output_bam: str, reference: str) -> None:
    """
    Anonymize the BAM file using GATK's ReadAnonymizer.

    Args:
        input_bam: Path to the input BAM file (after replacing read groups).
        output_bam: Path for the output anonymized BAM.
        reference: Path to the reference FASTA file.
    """
    logging.info("Running GATK ReadAnonymizer on: %s", input_bam)
    gatk_cmd = [
        "gatk", "ReadAnonymizer",
        "-I", input_bam,
        "-O", output_bam,
        "-R", reference,
    ]
    logging.debug("Running command: %s", " ".join(gatk_cmd))
    subprocess.run(gatk_cmd, check=True)


def run_samtools_reheader(input_bam: str, output_bam: str) -> None:
    """
    Reheader the BAM file to remove @PG, @RG, and @CO lines and index it.

    Args:
        input_bam: Path to the input BAM file (e.g. the anonymized BAM).
        output_bam: Path for the reheadered output BAM.
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


def remove_file_and_index(file_path: str) -> None:
    """
    Remove a file and any index files that may exist.
    Checks for both the standard naming scheme (file.bam.bai) and an alternative
    where the .bam extension is replaced by .bai.
    
    Args:
        file_path: Path to the file to remove.
    """
    try:
        if os.path.exists(file_path):
            os.remove(file_path)
            logging.debug("Removed file: %s", file_path)
    except Exception as e:
        logging.warning("Could not remove file %s: %s", file_path, e)

    # Check for possible index names.
    possible_indices = [file_path + ".bai", file_path.replace(".bam", ".bai")]
    for idx in possible_indices:
        try:
            if os.path.exists(idx):
                os.remove(idx)
                logging.debug("Removed index file: %s", idx)
        except Exception as e:
            logging.warning("Could not remove index file %s: %s", idx, e)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Subset a BAM file to a region, replace read groups with generic values, "
            "anonymize it using GATK ReadAnonymizer, and reheader to remove unwanted header lines. "
            "The final output file is named based on the provided target design."
        )
    )
    parser.add_argument("--input-bam", required=True,
                        help="Path to the input BAM file.")
    parser.add_argument("--target-design", required=True,
                        help="Target design name to be used for naming the output file (e.g. twist_v2).")
    parser.add_argument("--ref", required=True,
                        help="Path to the reference FASTA file (required by GATK tools).")
    parser.add_argument("--region", default="chr1:155184000-155194000",
                        help="Genomic region to subset from the BAM file. Default is 'chr1:155184000-155194000' (hg38 MUC1 VNTR region).")
    parser.add_argument("--output-dir", default=".",
                        help="Directory to write output files. Default is current directory.")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="Keep intermediate files for debugging purposes.")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Set logging level. Default is INFO.")
    return parser.parse_args()


def main() -> None:
    """Main function for the BAM anonymization pipeline."""
    args = parse_args()
    logging.basicConfig(level=args.log_level.upper(),
                        format="%(asctime)s [%(levelname)s] %(message)s")

    if not os.path.isfile(args.input_bam):
        logging.error("Input BAM file not found: %s", args.input_bam)
        sys.exit(1)
    if not os.path.isfile(args.ref):
        logging.error("Reference FASTA file not found: %s", args.ref)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Define intermediate and final file paths.
    subset_bam_path = os.path.join(args.output_dir, f"{args.target_design}_subset.bam")
    rg_bam_path = os.path.join(args.output_dir, f"{args.target_design}_rg.bam")
    anon_bam_path = os.path.join(args.output_dir, f"{args.target_design}_anon.bam")
    final_bam_path = os.path.join(args.output_dir, f"{args.target_design}.bam")

    try:
        # 1. Subset the input BAM to the specified region.
        subset_bam(args.input_bam, args.region, subset_bam_path)

        # 2. Replace read groups with generic values.
        run_gatk_add_or_replace_read_groups(subset_bam_path, rg_bam_path)

        # 3. Anonymize the BAM file.
        run_gatk_read_anonymizer(rg_bam_path, anon_bam_path, args.ref)

        # 4. Reheader the anonymized BAM.
        run_samtools_reheader(anon_bam_path, final_bam_path)

        logging.info("Final anonymized BAM file created successfully: %s", final_bam_path)

        # 5. Compute MD5 checksums and write pseudonymisation output.
        old_name = os.path.basename(args.input_bam)
        new_name = os.path.basename(final_bam_path)
        old_md5 = compute_md5(args.input_bam)
        new_md5 = compute_md5(final_bam_path)
        # Prefix the CSV file with the target design to avoid collisions.
        ps_output = os.path.join(args.output_dir, f"{args.target_design}_pseudonymisation_output.csv")
        with open(ps_output, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["old_name", "old_md5", "new_name", "new_md5"])
            writer.writerow([old_name, old_md5, new_name, new_md5])
        logging.info("Pseudonymisation output written to: %s", ps_output)
    except subprocess.CalledProcessError as err:
        logging.error("A subprocess failed: %s", err)
        sys.exit(1)
    finally:
        if args.keep_intermediates:
            logging.info("Keeping intermediate files (--keep-intermediates set).")
        else:
            for f in [subset_bam_path, rg_bam_path, anon_bam_path]:
                remove_file_and_index(f)


if __name__ == "__main__":
    main()

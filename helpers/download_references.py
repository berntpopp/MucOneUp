#!/usr/bin/env python3
"""
download_references.py

A helper script to download MUC1 reference sequences for both hg19 and hg38 assemblies.
This script:
1. Downloads left and right flanking regions of the MUC1 VNTR for both assemblies
2. Updates the config.json file with the downloaded sequences
3. Handles reverse complementing when necessary (MUC1 is on the negative strand)

Usage:
    python download_references.py [--config CONFIG_PATH]

Args:
    --config: Path to the configuration file (default: ../config.json)
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
import requests
import xml.etree.ElementTree as ET

# Set up basic logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

# Default padding around VNTR region (in base pairs)
DEFAULT_PADDING = 5000


def parse_region(region_str: str) -> tuple:
    """
    Parse a region string in the format "chrX:start-end" into components.

    Args:
        region_str: Region string in the format "chrX:start-end"

    Returns:
        tuple: (chromosome, start, end)
    """
    try:
        chrom, pos = region_str.split(":")
        start, end = map(int, pos.split("-"))
        return chrom, start, end
    except ValueError:
        logging.error(
            f"Invalid region format: {region_str}, expected format: chrX:start-end"
        )
        sys.exit(1)


def get_ucsc_region_das(
    assembly: str,
    chrom: str,
    start: int,  # DAS uses 1-based, inclusive coordinates
    end: int,  # DAS uses 1-based, inclusive coordinates
) -> str | None:
    """
    Fetches a specific genomic region from the UCSC DAS server.

    Args:
        assembly (str): Genome assembly name (e.g., 'hg38', 'hg19').
        chrom (str): Chromosome name (e.g., 'chr1', 'chrX').
        start (int): Start position (1-based, inclusive).
        end (int): End position (1-based, inclusive).

    Returns:
        str | None: The DNA sequence string if successful, None otherwise.
    """
    segment = f"{chrom}:{start},{end}"
    # Construct the DAS URL (adjust if needed for specific servers/assemblies)
    das_url = f"http://genome.ucsc.edu/cgi-bin/das/{assembly}/dna?segment={segment}"
    logging.info(f"Fetching sequence from: {das_url}")

    try:
        # Use a reasonable timeout for the request
        response = requests.get(das_url, timeout=30)
        response.raise_for_status()  # Check for HTTP errors (4xx or 5xx)

        # Parse the XML response
        root = ET.fromstring(response.content)

        # Find the DNA sequence within the XML structure
        # Structure is typically <DASDNA><SEQUENCE><DNA>...</DNA></SEQUENCE></DASDNA>
        dna_element = root.find(".//DNA")

        if dna_element is None or dna_element.text is None:
            logging.error(f"Could not find DNA sequence in DAS response for {segment}.")
            logging.debug(f"Response content:\n{response.text}")
            return None

        # The sequence might have newlines or other whitespace, remove it
        sequence = "".join(dna_element.text.split())
        # Convert sequence to uppercase
        sequence = sequence.upper()
        logging.info(
            f"Successfully fetched sequence for {segment} ({len(sequence)} bp)."
        )
        return sequence

    except requests.exceptions.Timeout:
        logging.error(f"Request timed out for {das_url}")
        return None
    except requests.exceptions.RequestException as e:
        logging.error(f"Error fetching data from UCSC DAS server: {e}")
        # Log details from the response if available
        if hasattr(e, "response") and e.response is not None:
            logging.error(f"Status Code: {e.response.status_code}")
            logging.error(f"Response Text: {e.response.text}")
        return None
    except ET.ParseError as e:
        logging.error(f"Error parsing XML response from DAS server: {e}")
        # Log the raw response text to help diagnose parsing issues
        if "response" in locals():
            logging.debug(f"Response content:\n{response.text}")
        return None
    except Exception as e:
        # Catch any other unexpected errors
        logging.error(
            f"An unexpected error occurred while fetching {segment}: {e}", exc_info=True
        )
        return None


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence

    Returns:
        str: Reverse complement of the sequence
    """
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "N": "N",
        "n": "n",
    }
    return "".join(complement.get(base, base) for base in reversed(sequence))


def download_flanking_regions(
    assembly: str, vntr_region: str, padding: int = DEFAULT_PADDING
) -> tuple:
    """
    Download the left and right flanking regions for a given VNTR region.

    Args:
        assembly: Genome assembly (hg19 or hg38)
        vntr_region: VNTR region in format "chrX:start-end"
        padding: Number of base pairs to include on each side

    Returns:
        tuple: (left_sequence, right_sequence)
    """
    chrom, start, end = parse_region(vntr_region)

    # MUC1 is on the negative strand, so we need to:
    # 1. Download the sequences
    # 2. Reverse complement them
    # 3. The "left" is actually to the right of the VNTR in genomic coordinates
    # 4. The "right" is actually to the left of the VNTR in genomic coordinates

    # Download left flanking region (right of VNTR in genomic coordinates)
    left_start = end + 1
    left_end = end + padding
    logging.info(
        f"Downloading {assembly} LEFT flanking region: {chrom}:{left_start}-{left_end} ({padding} bp)"
    )
    left_seq = get_ucsc_region_das(assembly, chrom, left_start, left_end)

    # Download right flanking region (left of VNTR in genomic coordinates)
    right_start = max(1, start - padding)
    right_end = start - 1
    logging.info(
        f"Downloading {assembly} RIGHT flanking region: {chrom}:{right_start}-{right_end} ({padding} bp)"
    )
    right_seq = get_ucsc_region_das(assembly, chrom, right_start, right_end)

    # Since MUC1 is on the negative strand, reverse complement both sequences
    if left_seq:
        logging.info(f"Reverse complementing LEFT sequence ({len(left_seq)} bp)")
        left_seq = reverse_complement(left_seq)
    if right_seq:
        logging.info(f"Reverse complementing RIGHT sequence ({len(right_seq)} bp)")
        right_seq = reverse_complement(right_seq)

    return left_seq, right_seq


def update_config_with_sequences(
    config_path: str, assembly: str, left_seq: str, right_seq: str
) -> bool:
    """
    Update the config file with the downloaded sequences.

    Args:
        config_path: Path to the config JSON file
        assembly: Assembly name (hg19 or hg38)
        left_seq: Left flanking sequence
        right_seq: Right flanking sequence

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with open(config_path, "r") as f:
            config = json.load(f)

        # Update the sequences in the constants section
        if assembly in config["constants"]:
            config["constants"][assembly]["left"] = left_seq
            config["constants"][assembly]["right"] = right_seq

            # Write back the updated config
            with open(config_path, "w") as f:
                json.dump(config, f, indent=2)

            logging.info(f"Updated config with {assembly} sequences")
            return True
        else:
            logging.error(f"Assembly {assembly} not found in config constants section")
            return False

    except Exception as e:
        logging.error(f"Error updating config: {e}")
        return False


def main():
    """Main function to download references and update config."""
    parser = argparse.ArgumentParser(
        description="Download reference sequences for MUC1 region in hg19 and hg38."
    )
    parser.add_argument(
        "--config",
        default=os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config.json"
        ),
        help="Path to the configuration file",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=DEFAULT_PADDING,
        help=f"Padding around VNTR region in bp (default: {DEFAULT_PADDING})",
    )
    parser.add_argument(
        "--assembly",
        choices=["hg19", "hg38", "both"],
        default="both",
        help="Which assembly to download references for (default: both)",
    )

    args = parser.parse_args()

    # Check if the config file exists
    if not os.path.exists(args.config):
        logging.error(f"Config file not found: {args.config}")
        sys.exit(1)

    # Load the config to get the VNTR regions
    with open(args.config, "r") as f:
        config = json.load(f)

    assemblies = []
    if args.assembly == "both":
        assemblies = ["hg19", "hg38"]
    else:
        assemblies = [args.assembly]

    for assembly in assemblies:
        try:
            vntr_region_key = f"vntr_region_{assembly}"
            vntr_region = config["read_simulation"][vntr_region_key]

            logging.info(
                f"Downloading {assembly} flanking regions for VNTR: {vntr_region}"
            )
            left_seq, right_seq = download_flanking_regions(
                assembly, vntr_region, args.padding
            )

            if left_seq and right_seq:
                logging.info(
                    f"Downloaded {assembly} sequences: LEFT={len(left_seq)}bp, RIGHT={len(right_seq)}bp"
                )
                update_config_with_sequences(args.config, assembly, left_seq, right_seq)
            else:
                logging.error(f"Failed to download {assembly} sequences")
        except KeyError as e:
            logging.error(f"Missing key in config: {e}")
        except Exception as e:
            logging.error(f"Error processing {assembly}: {e}")

    logging.info("Reference download process completed.")


if __name__ == "__main__":
    # Ensure required libraries are available
    try:
        import requests
        import xml.etree.ElementTree
    except ImportError as e:
        print(f"Error: Required library not found: {e.name}", file=sys.stderr)
        print("Please install it using: pip install requests", file=sys.stderr)
        sys.exit(1)

    main()

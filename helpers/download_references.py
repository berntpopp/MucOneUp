#!/usr/bin/env python3
"""
download_references.py

A helper script to download MUC1 reference sequences for both hg19 and hg38 assemblies,
and set up the required references for both Illumina and Oxford Nanopore sequencing.

This script:
1. Downloads left and right flanking regions of the MUC1 VNTR for both assemblies
2. Updates the config.json file with the downloaded sequences
3. Handles reverse complementing when necessary (MUC1 is on the negative strand)
4. Sets up NanoSim model for Oxford Nanopore simulation
5. Creates minimap2 indices for human reference genomes

Usage:
    python download_references.py [--config CONFIG_PATH]

Args:
    --config: Path to the configuration file (default: ../config.json)
    --references-dir: Path to the references directory (default: ../reference)
    --setup-nanosim: Flag to download and set up NanoSim models
    --create-minimap2-index: Flag to create minimap2 indices for human references
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tarfile
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


def download_nanosim_model(output_dir, force=False):
    """Download and extract the NanoSim model.

    Args:
        output_dir: Directory to store the NanoSim model
        force: Force download even if model exists

    Returns:
        Path to the extracted model directory
    """
    # NanoSim model URL
    model_url = "https://github.com/bcgsc/NanoSim/raw/master/pre-trained_models/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1.tar.gz"
    model_name = "human_giab_hg002_sub1M_kitv14_dorado_v3.2.1"

    # Create nanosim directory in output_dir
    nanosim_dir = Path(output_dir) / "nanosim"
    nanosim_dir.mkdir(parents=True, exist_ok=True)

    # Paths for downloaded and extracted files
    tar_path = nanosim_dir / f"{model_name}.tar.gz"
    extract_path = nanosim_dir / model_name

    # Check if model already exists
    if extract_path.exists() and not force:
        logging.info(f"NanoSim model already exists at {extract_path}")
        return extract_path

    # Download the model
    logging.info(f"Downloading NanoSim model from {model_url}")
    try:
        response = requests.get(model_url, stream=True)
        response.raise_for_status()
        with open(tar_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        logging.info(f"Downloaded NanoSim model to {tar_path}")
    except Exception as e:
        logging.error(f"Failed to download NanoSim model: {e}")
        return None

    # Extract the model
    try:
        if extract_path.exists() and force:
            logging.info(f"Removing existing model directory: {extract_path}")
            shutil.rmtree(extract_path)

        extract_path.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path) as tar:
            tar.extractall(path=extract_path)
        logging.info(f"Extracted NanoSim model to {extract_path}")
        return extract_path
    except Exception as e:
        logging.error(f"Failed to extract NanoSim model: {e}")
        return None


def create_minimap2_index(reference_path, output_dir, minimap2_cmd="minimap2"):
    """Create a minimap2 index for the reference genome.

    Args:
        reference_path: Path to the reference genome FASTA
        output_dir: Directory to store the index
        minimap2_cmd: Command to run minimap2 (default: 'minimap2')

    Returns:
        Path to the created index file
    """
    # Check if reference exists
    ref_path = Path(reference_path)
    if not ref_path.exists():
        logging.error(f"Reference file not found: {reference_path}")
        return None

    # Create index path
    index_path = Path(output_dir) / f"{ref_path.name}.mmi"

    # Check if index already exists
    if index_path.exists():
        logging.info(f"Minimap2 index already exists at {index_path}")
        return index_path

    # Create the index
    logging.info(f"Creating minimap2 index for {reference_path}")
    try:
        # Convert paths to strings to avoid any path parsing issues
        ref_path_str = str(ref_path)
        index_path_str = str(index_path)

        # Construct the command
        # We use "mamba run" with env_nanosim to ensure correct environment
        cmd = f"mamba run --no-capture-output -n env_nanosim {minimap2_cmd} -d {index_path_str} {ref_path_str}"
        logging.info(f"Running command: {cmd}")

        # Run the command
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        logging.info(f"Successfully created minimap2 index at {index_path}")
        return index_path
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create minimap2 index: {e}")
        logging.error(f"Command output: {e.stdout}")
        logging.error(f"Command error: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Error creating minimap2 index: {e}")
        return None


def main():
    """Main function to download references and update config."""
    parser = argparse.ArgumentParser(
        description="Download reference sequences for MUC1 region in hg19 and hg38 and set up NanoSim models."
    )
    parser.add_argument(
        "--config",
        default=os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config.json"
        ),
        help="Path to the configuration file",
    )
    parser.add_argument(
        "--references-dir",
        default=os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference"
        ),
        help="Path to the references directory",
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
    parser.add_argument(
        "--download-flanking",
        action="store_true",
        help="Download flanking regions for MUC1 VNTR",
    )
    parser.add_argument(
        "--setup-nanosim",
        action="store_true",
        help="Download and set up NanoSim models",
    )
    parser.add_argument(
        "--create-minimap2-index",
        action="store_true",
        help="Create minimap2 indices for human references",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force download and extraction even if files exist",
    )

    args = parser.parse_args()

    # Check if the config file exists
    if not os.path.exists(args.config):
        logging.error(f"Config file not found: {args.config}")
        sys.exit(1)

    # Check if any specific operations are requested
    specific_ops_requested = (
        args.download_flanking or args.setup_nanosim or args.create_minimap2_index
    )

    # Load the config
    with open(args.config, "r") as f:
        config = json.load(f)

    # Process assembly choice
    assemblies = []
    if args.assembly == "both":
        assemblies = ["hg19", "hg38"]
    else:
        assemblies = [args.assembly]

    # Download flanking regions if explicitly requested or if no specific operations are requested
    if args.download_flanking or not specific_ops_requested:
        for assembly in assemblies:
            try:
                # Try to get the vntr_region from constants first
                if "vntr_region" in config["constants"][assembly]:
                    vntr_region = config["constants"][assembly]["vntr_region"]
                else:
                    # Fall back to read_simulation section
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
                    update_config_with_sequences(
                        args.config, assembly, left_seq, right_seq
                    )
                else:
                    logging.error(f"Failed to download {assembly} sequences")
            except KeyError as e:
                logging.error(f"Missing key in config: {e}")
            except Exception as e:
                logging.error(f"Error processing {assembly}: {e}")

    # Setup NanoSim models if requested
    if args.setup_nanosim or not specific_ops_requested:
        logging.info("Setting up NanoSim models...")
        nanosim_model_path = download_nanosim_model(args.references_dir, args.force)
        if nanosim_model_path:
            # Update the config with the NanoSim model path
            try:
                with open(args.config, "r") as f:
                    config = json.load(f)

                # Update the NanoSim training model path in config
                if "nanosim_params" in config:
                    rel_path = os.path.relpath(
                        nanosim_model_path, os.path.dirname(args.config)
                    )
                    config["nanosim_params"]["training_data_path"] = str(rel_path)
                    logging.info(f"Updated config with NanoSim model path: {rel_path}")

                    # Write the updated config
                    with open(args.config, "w") as f:
                        json.dump(config, f, indent=2)
                else:
                    logging.warning(
                        "nanosim_params section not found in config, could not update model path"
                    )
            except Exception as e:
                logging.error(f"Error updating config with NanoSim model path: {e}")

    # Create minimap2 indices if requested
    if args.create_minimap2_index or not specific_ops_requested:
        logging.info("Creating minimap2 indices for human references...")
        try:
            with open(args.config, "r") as f:
                config = json.load(f)

            # Get the human reference path from config
            human_reference = config.get("read_simulation", {}).get("human_reference")
            if human_reference:
                # Convert to absolute path
                if not os.path.isabs(human_reference):
                    human_reference = os.path.join(
                        os.path.dirname(args.config), human_reference
                    )

                # Ensure the path exists
                if not os.path.exists(human_reference):
                    logging.error(f"Human reference file not found: {human_reference}")
                else:
                    # Create the index
                    index_path = create_minimap2_index(
                        human_reference, os.path.dirname(human_reference)
                    )
                if index_path:
                    logging.info(f"Successfully created minimap2 index at {index_path}")
            else:
                logging.warning(
                    "human_reference not found in config, could not create minimap2 index"
                )
        except Exception as e:
            logging.error(f"Error creating minimap2 index: {e}")

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

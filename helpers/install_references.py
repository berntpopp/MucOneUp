#!/usr/bin/env python3
"""
Install Reference Sequences for muconeup

This script downloads and installs external reference files required for the muconeup
pipeline. It downloads files specified in a configuration file, verifies file integrity
via MD5, extracts compressed files if needed, optionally indexes FASTA files using BWA,
and (if configured) generates a sequence dictionary using GATK.

Usage:
    python install_references.py --output-dir /path/to/destination [--config /path/to/config.json] [--skip-indexing]

Example:
    python install_references.py --output-dir ./references

The default configuration file is 'install_references_config.json' (located alongside this script).
A sample configuration is provided below.
"""

import argparse
import gzip
import hashlib
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

from urllib.request import urlretrieve


def report_progress(block_num: int, block_size: int, total_size: int) -> None:
    """
    Report download progress to stdout.

    Args:
        block_num (int): Number of blocks transferred so far.
        block_size (int): Size of each block (in bytes).
        total_size (int): Total size of the file (in bytes).
    """
    downloaded = block_num * block_size
    if total_size > 0:
        percent = downloaded / total_size * 100
        if percent > 100:
            percent = 100
        sys.stdout.write(f"\rDownloading... {downloaded} of {total_size} bytes ({percent:.1f}%)")
    else:
        sys.stdout.write(f"\rDownloading... {downloaded} bytes")
    sys.stdout.flush()


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load JSON configuration from a file.

    Args:
        config_path (Path): Path to the configuration JSON file.

    Returns:
        Dict[str, Any]: Configuration dictionary.

    Exits:
        Exits with error if file not found or JSON is invalid.
    """
    if not config_path.exists():
        logging.error("Configuration file not found at %s", config_path)
        sys.exit(1)
    try:
        with config_path.open("r") as file:
            config = json.load(file)
        return config
    except json.JSONDecodeError as err:
        logging.error("Error parsing JSON configuration: %s", err)
        sys.exit(1)
    except Exception as err:
        logging.error("Unexpected error loading configuration: %s", err)
        sys.exit(1)


def download_file(url: str, dest_path: Path) -> None:
    """Download a file from a URL if it does not already exist.

    Args:
        url (str): URL to download from.
        dest_path (Path): Destination file path.
    """
    if dest_path.exists():
        logging.info("File already exists at %s. Skipping download.", dest_path)
        return

    logging.info("Downloading %s to %s", url, dest_path)
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        # Show progress via report_progress
        urlretrieve(url, dest_path, reporthook=report_progress)
        # Ensure a newline after download completes
        sys.stdout.write("\n")
        logging.info("Downloaded file: %s", dest_path.name)
    except Exception as err:
        logging.error("Failed to download %s: %s", url, err)
        sys.exit(1)


def calculate_md5(file_path: Path, chunk_size: int = 4096) -> str:
    """Calculate the MD5 checksum of a file.

    Args:
        file_path (Path): File to checksum.
        chunk_size (int): Bytes to read per iteration.

    Returns:
        str: MD5 checksum.
    """
    md5_hash = hashlib.md5()
    try:
        with file_path.open("rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                md5_hash.update(chunk)
        checksum = md5_hash.hexdigest()
        logging.info("MD5 checksum for %s: %s", file_path.name, checksum)
        return checksum
    except Exception as err:
        logging.error("Error calculating MD5 for %s: %s", file_path, err)
        sys.exit(1)


def extract_gzip(file_path: Path) -> Path:
    """Extract a gzip-compressed file.

    Args:
        file_path (Path): Path to the .gz file.

    Returns:
        Path: Path to the extracted file.
    """
    extracted_path = file_path.with_suffix("")
    logging.info("Extracting %s to %s", file_path.name, extracted_path.name)
    try:
        with gzip.open(file_path, "rb") as f_in, extracted_path.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logging.info("Extraction complete: %s", extracted_path.name)
        return extracted_path
    except Exception as err:
        logging.error("Failed to extract %s: %s", file_path, err)
        sys.exit(1)


def execute_index_command(command_template: str, fasta_path: Path) -> None:
    """Execute a command on a FASTA file.

    Args:
        command_template (str): Command template with a {path} placeholder.
        fasta_path (Path): Path to the FASTA file.
    """
    command = command_template.format(path=str(fasta_path))
    logging.info("Executing command: %s", command)
    try:
        subprocess.run(
            command.split(),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        logging.info("Command completed for %s", fasta_path.name)
    except subprocess.CalledProcessError as err:
        logging.error("Command failed for %s: %s", fasta_path, err.stderr.decode().strip())
        sys.exit(1)
    except Exception as err:
        logging.error("Error executing command for %s: %s", fasta_path, err)
        sys.exit(1)


def setup_logging(output_dir: Path) -> None:
    """Set up logging to stdout and to a log file in the output directory.

    Args:
        output_dir (Path): Directory to place the log file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "install_references.log"

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Remove existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    stream_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(log_file)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    stream_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logging.info("Logging initialized. Log file: %s", log_file)


def write_installed_config(installed_refs: Dict[str, str], output_dir: Path) -> None:
    """Write a JSON file mapping reference names to installed file paths.

    Args:
        installed_refs (Dict[str, str]): Mapping of reference names to file paths.
        output_dir (Path): Directory to save the configuration file.
    """
    output_file = output_dir / "installed_references.json"
    try:
        with output_file.open("w") as f:
            json.dump(installed_refs, f, indent=2)
        logging.info("Installed references configuration written to %s", output_file)
    except Exception as err:
        logging.error("Failed to write installed references configuration: %s", err)
        sys.exit(1)


def process_reference(
    ref_name: str,
    ref_info: Dict[str, Any],
    output_dir: Path,
    skip_indexing: bool,
    bwa_path: str
) -> str:
    """Download, verify, extract, and index (if needed) a single reference.

    Args:
        ref_name (str): Identifier for the reference.
        ref_info (Dict[str, Any]): Reference details (url, target_path, index_command, seq_dict_command).
        output_dir (Path): Base directory for downloads.
        skip_indexing (bool): Flag to skip indexing.
        bwa_path (str): Path to the BWA executable.

    Returns:
        str: Absolute path to the installed reference file.
    """
    url = ref_info.get("url")
    target_rel_path = ref_info.get("target_path")
    if not url or not target_rel_path:
        logging.error("Missing 'url' or 'target_path' for reference %s", ref_name)
        sys.exit(1)

    target_path = output_dir / target_rel_path
    download_file(url, target_path)
    calculate_md5(target_path)

    # If the file is gzip-compressed, extract it.
    if target_path.suffix == ".gz":
        installed_path = extract_gzip(target_path)
    else:
        installed_path = target_path

    # Run indexing command if provided and not skipped.
    index_command = ref_info.get("index_command")
    if index_command and not skip_indexing:
        # Replace "bwa" with the executable from the config if needed.
        command = index_command.replace("bwa", bwa_path)
        execute_index_command(command, installed_path)
    elif index_command:
        logging.info("Skipping indexing for %s", ref_name)

    # Run sequence dictionary command if provided.
    seq_dict_command = ref_info.get("seq_dict_command")
    if seq_dict_command:
        logging.info("Executing sequence dictionary command for %s", ref_name)
        execute_index_command(seq_dict_command, installed_path)

    return str(installed_path.resolve())


def main() -> None:
    """Main routine to download and install references."""
    parser = argparse.ArgumentParser(
        description="Download and install reference sequences for muconeup."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Destination directory for installed reference files."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "install_references_config.json",
        help="Path to the configuration JSON file."
    )
    parser.add_argument(
        "--skip-indexing",
        action="store_true",
        help="Skip indexing of FASTA reference files."
    )
    args = parser.parse_args()

    setup_logging(args.output_dir)
    config = load_config(args.config)
    bwa_path = config.get("bwa_path", "bwa")
    references = config.get("references", {})

    if not references:
        logging.error("No references found in configuration.")
        sys.exit(1)

    installed_refs = {}
    for ref_name, ref_info in references.items():
        logging.info("Processing reference: %s", ref_name)
        installed_path = process_reference(
            ref_name, ref_info, args.output_dir, args.skip_indexing, bwa_path
        )
        installed_refs[ref_name] = installed_path

    write_installed_config(installed_refs, args.output_dir)
    logging.info("All references have been installed successfully.")


if __name__ == "__main__":
    main()

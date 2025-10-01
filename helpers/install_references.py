#!/usr/bin/env python3
"""
Install Reference Sequences for muconeup

This script downloads and installs external reference files required for the muconeup
pipeline. It downloads files specified in a configuration file, verifies file integrity
via MD5, extracts compressed files if needed, optionally indexes FASTA files using BWA,
and (if configured) generates a sequence dictionary using GATK.

Usage:
    python install_references.py --output-dir /path/to/destination [--config /path/to/config.json] [--skip-indexing] [--force]

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
import tarfile
from pathlib import Path
from typing import Any
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


def load_config(config_path: Path) -> dict[str, Any]:
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
    """Download a file from a URL.

    Args:
        url (str): URL to download from.
        dest_path (Path): Destination file path.
    """
    logging.info("Downloading %s to %s", url, dest_path)
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(url, dest_path, reporthook=report_progress)
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
        subprocess.run(command.split(), check=True, capture_output=True)
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


def write_installed_config(installed_refs: dict[str, str], output_dir: Path) -> None:
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
    ref_info: dict[str, Any],
    output_dir: Path,
    skip_indexing: bool,
    bwa_path: str,
    gatk_path: str,
    minimap2_path: str,
    force: bool,
) -> str:
    """Download, verify, extract, index (if needed), and generate a sequence dictionary (if configured) for a single reference.

    Args:
        ref_name (str): Identifier for the reference.
        ref_info (Dict[str, Any]): Reference details (url, target_path, index_command, seq_dict_command).
        output_dir (Path): Base directory for downloads.
        skip_indexing (bool): Flag to skip indexing.
        bwa_path (str): Path to the BWA executable.
        gatk_path (str): Path to the GATK executable.
        force (bool): If True, re-run all steps regardless of existing files.

    Returns:
        str: Absolute path to the installed reference file.
    """
    url = ref_info.get("url")
    target_rel_path = ref_info.get("target_path")
    if not url or not target_rel_path:
        logging.error("Missing 'url' or 'target_path' for reference %s", ref_name)
        sys.exit(1)

    target_path = output_dir / target_rel_path
    if force and target_path.exists():
        logging.info("Force enabled, removing existing file: %s", target_path)
        target_path.unlink()

    download_file(url, target_path)
    calculate_md5(target_path)

    # Check if this is a tar file that needs extraction
    if ref_info.get("extract_tar", False) and (
        str(target_path).endswith(".tar.gz") or str(target_path).endswith(".tar")
    ):
        extract_path = ref_info.get("extract_path")
        if not extract_path:
            extract_path = target_path.stem
            if extract_path.endswith(".tar"):
                extract_path = extract_path[:-4]  # Remove .tar suffix if present

        extract_dir = output_dir / extract_path
        if extract_dir.exists() and not force:
            logging.info(
                "Extracted directory %s already exists, skipping extraction.",
                extract_dir,
            )
            installed_path = extract_dir
        else:
            if force and extract_dir.exists() and extract_dir.is_dir():
                logging.info(
                    "Force enabled, removing existing extracted directory: %s",
                    extract_dir,
                )
                shutil.rmtree(extract_dir)
            # Extract the tar file
            with tarfile.open(target_path) as tar:
                extract_dir.mkdir(parents=True, exist_ok=True)
                logging.info("Extracting %s to %s", target_path, extract_dir)
                tar.extractall(path=extract_dir)
            installed_path = extract_dir
            logging.info("Successfully extracted %s to %s", target_path.name, extract_dir)
    # If the file is gzip-compressed, check for extracted file and extract if needed.
    elif target_path.suffix == ".gz":
        extracted_path = target_path.with_suffix("")
        if extracted_path.exists() and not force:
            logging.info(
                "Extracted file %s already exists, skipping extraction.",
                extracted_path.name,
            )
            installed_path = extracted_path
        else:
            if force and extracted_path.exists():
                logging.info(
                    "Force enabled, removing existing extracted file: %s",
                    extracted_path,
                )
                extracted_path.unlink()
            installed_path = extract_gzip(target_path)
    else:
        installed_path = target_path

    # Run indexing command if provided and not skipped.
    index_command = ref_info.get("index_command")
    if index_command and not skip_indexing:
        # Determine one of the expected index files (e.g. .amb).
        index_file = installed_path.parent / (installed_path.name + ".amb")
        if index_file.exists() and not force:
            logging.info(
                "Index files already exist for %s. Skipping indexing.",
                installed_path.name,
            )
        else:
            if force and index_file.exists():
                # Remove all expected index files.
                for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
                    f = installed_path.parent / (installed_path.name + ext)
                    if f.exists():
                        logging.info("Force enabled, removing existing index file: %s", f)
                        f.unlink()
            command = index_command.replace("bwa", bwa_path)
            execute_index_command(command, installed_path)
    elif index_command:
        logging.info("Skipping indexing for %s", ref_name)

    # Run minimap2 indexing command if provided and not skipped
    minimap2_index_command = ref_info.get("minimap2_index_command")
    if minimap2_index_command and not skip_indexing:
        # Check for .mmi index file
        index_file = Path(str(installed_path) + ".mmi")
        if index_file.exists() and not force:
            logging.info(
                "Minimap2 index file already exists for %s. Skipping indexing.",
                installed_path.name,
            )
        else:
            if force and index_file.exists():
                logging.info(
                    "Force enabled, removing existing minimap2 index file: %s",
                    index_file,
                )
                index_file.unlink()
            command = minimap2_index_command.replace("minimap2", minimap2_path)
            logging.info("Executing minimap2 indexing command for %s", ref_name)
            execute_index_command(command, installed_path)
    elif minimap2_index_command:
        logging.info("Skipping minimap2 indexing for %s", ref_name)

    # Run sequence dictionary command if provided.
    seq_dict_command = ref_info.get("seq_dict_command")
    if seq_dict_command:
        seq_dict_file = installed_path.parent / (installed_path.name + ".dict")
        if seq_dict_file.exists() and not force:
            logging.info(
                "Sequence dictionary already exists for %s. Skipping creation.",
                installed_path.name,
            )
        else:
            if force and seq_dict_file.exists():
                logging.info(
                    "Force enabled, removing existing sequence dictionary: %s",
                    seq_dict_file,
                )
                seq_dict_file.unlink()
            logging.info("Executing sequence dictionary command for %s", ref_name)
            seq_dict_command = seq_dict_command.replace("gatk", gatk_path)
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
        help="Destination directory for installed reference files.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "install_references_config.json",
        help="Path to the configuration JSON file.",
    )
    parser.add_argument(
        "--skip-indexing",
        action="store_true",
        help="Skip indexing of FASTA reference files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Redo everything despite existing files and overwrite them.",
    )
    args = parser.parse_args()

    setup_logging(args.output_dir)
    config = load_config(args.config)
    bwa_path = config.get("bwa_path", "bwa")
    gatk_path = config.get("gatk_path", "gatk")
    minimap2_path = config.get("minimap2_path", "minimap2")
    references = config.get("references", {})

    if not references:
        logging.error("No references found in configuration.")
        sys.exit(1)

    installed_refs = {}
    for ref_name, ref_info in references.items():
        logging.info("Processing reference: %s", ref_name)
        installed_path = process_reference(
            ref_name,
            ref_info,
            args.output_dir,
            args.skip_indexing,
            bwa_path,
            gatk_path,
            minimap2_path,
            args.force,
        )
        installed_refs[ref_name] = installed_path

    write_installed_config(installed_refs, args.output_dir)
    logging.info("All references have been installed successfully.")


if __name__ == "__main__":
    main()

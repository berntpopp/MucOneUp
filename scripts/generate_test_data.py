#!/usr/bin/env python3
"""
Generate comprehensive MucOneUp test dataset.

This script generates test data across multiple sequencing platforms
with clean code architecture following SOLID principles.

Usage:
    python scripts/generate_test_data.py --version v0.25.0 --config config.json

For help:
    python scripts/generate_test_data.py --help
"""

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

# Add testdata module to path
sys.path.insert(0, str(Path(__file__).parent))

from testdata.config import load_config
from testdata.generator import ReadSimulator, ReferenceGenerator
from testdata.metadata import MetadataGenerator
from testdata.packaging import Packager
from testdata.utils import collect_all_files, ensure_dir


def setup_logging(verbose: bool = False):
    """
    Configure logging for the script.

    Args:
        verbose: Enable debug-level logging
    """
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("generate_test_data.log"),
        ],
    )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate MucOneUp test dataset with multi-platform support"
    )
    parser.add_argument("--version", default="v0.25.0", help="Version tag (default: v0.25.0)")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("test_data"),
        help="Output directory (default: test_data/)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.json"),
        help="MucOneUp config file (default: config.json)",
    )
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument(
        "--platforms",
        nargs="+",
        choices=["illumina", "ont", "pacbio", "all"],
        default=["all"],
        help="Platforms to generate (default: all)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose debug logging")
    parser.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Keep intermediate files (for debugging)",
    )
    parser.add_argument("--no-tarball", action="store_true", help="Skip tarball creation")

    return parser.parse_args()


def main():
    """Main entry point for test data generation."""
    args = parse_args()

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Validate config file exists
    if not args.config.exists():
        logger.error(f"Config file not found: {args.config}")
        return 1

    # Load configuration
    config = load_config(version=args.version, threads=args.threads, verbose=args.verbose)

    # Filter platforms if specified
    if "all" not in args.platforms:
        for name in list(config.platforms.keys()):
            if name not in args.platforms:
                config.platforms[name].enabled = False

    # Create output directory
    dataset_dir = args.output / f"testdata_40-70_{args.version}"
    ensure_dir(dataset_dir)

    logger.info("=" * 70)
    logger.info("MucOneUp Test Dataset Generator v2.0")
    logger.info("=" * 70)
    logger.info(f"Version: {args.version}")
    logger.info(f"Output: {dataset_dir}")
    logger.info(f"Platforms: {[p.name for p in config.enabled_platforms]}")
    logger.info(f"Threads: {args.threads}")
    logger.info("=" * 70)

    try:
        # Step 1: Generate diploid references
        logger.info("")
        logger.info("STEP 1: Generating diploid references")
        logger.info("-" * 70)
        ref_generator = ReferenceGenerator(config, args.config)
        references = ref_generator.generate(dataset_dir)

        # Step 2: Simulate reads for each platform
        logger.info("")
        logger.info("STEP 2: Simulating reads for enabled platforms")
        logger.info("-" * 70)
        read_simulator = ReadSimulator(config, args.config)

        platform_files = {}
        for platform in config.enabled_platforms:
            files = read_simulator.simulate_platform(platform, references, dataset_dir)
            platform_files[platform.name] = files

        # Step 3: Generate metadata
        logger.info("")
        logger.info("STEP 3: Generating metadata")
        logger.info("-" * 70)

        # Collect all files
        all_files = collect_all_files(
            dataset_dir, exclude_patterns=["*.log", "*.pyc", "__pycache__"]
        )

        metadata_gen = MetadataGenerator(config)

        manifest_path = metadata_gen.generate_manifest(dataset_dir, all_files)

        # Load manifest for checksums file
        with open(manifest_path) as f:
            manifest = json.load(f)

        metadata_gen.generate_checksums_file(dataset_dir, manifest)
        metadata_gen.generate_readme(dataset_dir)

        # Copy config for reproducibility
        shutil.copy(args.config, dataset_dir / "config.json")
        logger.info("✓ Copied config.json for reproducibility")

        # Step 4: Create tarball (optional)
        tarball_path = None
        if not args.no_tarball:
            logger.info("")
            logger.info("STEP 4: Creating tarball")
            logger.info("-" * 70)
            packager = Packager(args.version)
            tarball_path = packager.create_tarball(dataset_dir)

        # Cleanup if requested
        if not args.no_cleanup:
            logger.info("")
            logger.info("Cleaning up (use --no-cleanup to keep all files)")
            # Could add cleanup logic here if needed

        # Final summary
        logger.info("")
        logger.info("=" * 70)
        logger.info("✅ Test dataset generation complete!")
        logger.info("=" * 70)
        logger.info(f"📁 Dataset: {dataset_dir}")
        if tarball_path:
            logger.info(f"📦 Tarball: {tarball_path.name}")
            logger.info(f"   Size: {tarball_path.stat().st_size / 1024 / 1024:.2f} MB")
        logger.info(f"📋 Files: {len(all_files)}")
        logger.info(f"🧬 Platforms: {len(config.enabled_platforms)}")
        logger.info("🔐 Checksums: checksums.txt")
        logger.info("=" * 70)
        if tarball_path:
            logger.info(f"Next: Upload {tarball_path.name} to GitHub Release {args.version}")
        logger.info("=" * 70)

        return 0

    except Exception as e:
        logger.error(f"❌ Generation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

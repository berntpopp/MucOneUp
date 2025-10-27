#!/usr/bin/env python3
"""
Generate comprehensive MucOneUp test dataset.

This script generates test data across multiple sequencing platforms
(Illumina, ONT, PacBio) with unified logging, ORF analysis, and SNaPshot validation.

Usage:
    python helpers/generate_test_data.py --version v0.25.0 --config config.json

For help:
    python helpers/generate_test_data.py --help
"""

import argparse
import hashlib
import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, cast

# ============================================================================
# Configuration
# ============================================================================

DEFAULT_VERSION = "v0.25.0"
DEFAULT_OUTPUT = Path("helpers/test_data")  # Temporary working directory
DEFAULT_CONFIG = Path("config.json")
DEFAULT_THREADS = 4

# Sample configuration
SAMPLE_BASE_NAME = "testdata_40-70"
VNTR_H1_LENGTH = 40  # Haplotype 1 VNTR length
VNTR_H2_LENGTH = 70  # Haplotype 2 VNTR length
MUTATION_NAME = "normal,dupC"  # Dual simulation mode
MUTATION_TARGETS = "1,20"  # dupC at position 20 of haplotype 1
SEED_BASE = 42000

# Platform configurations
PLATFORMS = {
    "illumina": {"enabled": True, "coverage": 50, "seed_offset": 0, "output_subdir": "illumina"},
    "ont": {"enabled": True, "coverage": 30, "seed_offset": 10, "output_subdir": "ont"},
    "pacbio": {"enabled": True, "coverage": 30, "seed_offset": 20, "output_subdir": "pacbio"},
}

# ============================================================================
# Logging Setup
# ============================================================================


def setup_logging(verbose: bool = False):
    """Configure logging for the script."""
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("generate_test_data.log"),
        ],
    )


# ============================================================================
# Helper Functions
# ============================================================================


def ensure_dir(path: Path) -> Path:
    """Ensure directory exists, create if needed."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_command(
    cmd: list[str], description: str, check: bool = True, log_file: Path | None = None
) -> subprocess.CompletedProcess:
    """Run a command with logging."""
    logger = logging.getLogger(__name__)

    logger.info(f"{description}...")
    logger.debug(f"Command: {' '.join(str(c) for c in cmd)}")

    try:
        # If log_file specified, write command to it
        if log_file:
            with open(log_file, "a") as f:
                f.write(f"\n{'=' * 80}\n")
                f.write(f"Command: {' '.join(str(c) for c in cmd)}\n")
                f.write(f"Description: {description}\n")
                f.write(f"Timestamp: {datetime.now().isoformat()}\n")
                f.write(f"{'=' * 80}\n\n")

        result = subprocess.run(cmd, check=check, capture_output=True, text=True)

        # Write output to log file if specified
        if log_file:
            with open(log_file, "a") as f:
                if result.stdout:
                    f.write("STDOUT:\n")
                    f.write(result.stdout)
                    f.write("\n")
                if result.stderr:
                    f.write("STDERR:\n")
                    f.write(result.stderr)
                    f.write("\n")
                f.write(f"Exit code: {result.returncode}\n\n")

        if result.stdout:
            logger.debug(f"stdout: {result.stdout[:500]}")
        if result.stderr:
            logger.debug(f"stderr: {result.stderr[:500]}")

        logger.info(f"✓ {description} completed")

        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"✗ {description} failed")
        logger.error(f"Exit code: {e.returncode}")

        # Write error to log file
        if log_file:
            with open(log_file, "a") as f:
                f.write(f"ERROR: Command failed with exit code {e.returncode}\n")
                if e.stdout:
                    f.write("STDOUT:\n")
                    f.write(e.stdout)
                    f.write("\n")
                if e.stderr:
                    f.write("STDERR:\n")
                    f.write(e.stderr)
                    f.write("\n\n")

        if e.stdout:
            logger.error(f"stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"stderr: {e.stderr}")
        raise


def calculate_checksum(file_path: Path) -> str:
    """Calculate SHA256 checksum for a file."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def collect_all_files(root_dir: Path, exclude_patterns: list[str] | None = None) -> list[Path]:
    """Recursively collect all files in directory."""
    exclude_patterns = exclude_patterns or []
    files = []

    for item in root_dir.rglob("*"):
        if item.is_file():
            # Check exclusions
            skip = False
            for pattern in exclude_patterns:
                if item.match(pattern):
                    skip = True
                    break

            if not skip:
                files.append(item)

    return sorted(files)


# ============================================================================
# Tool Version Capture
# ============================================================================


def capture_tool_versions(output_dir: Path):
    """Capture versions of all tools."""
    logger = logging.getLogger(__name__)

    version_file = output_dir / "tool_versions.txt"
    with open(version_file, "w") as f:
        f.write(f"Tool Versions - Generated: {datetime.now().isoformat()}\n")
        f.write(f"{'=' * 70}\n\n")

        for tool, cmd in [
            ("MucOneUp", ["muconeup", "--version"]),
            ("Python", ["python", "--version"]),
            ("BWA", ["bwa"]),
            ("Samtools", ["samtools", "--version"]),
        ]:
            f.write(f"{tool}:\n")
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
                output = (result.stdout or result.stderr).strip().split("\n")[0]
                f.write(f"  {output}\n")
            except Exception:
                f.write("  unknown\n")
            f.write("\n")

    logger.info(f"✓ Tool versions captured: {version_file}")


# ============================================================================
# Reference Generation
# ============================================================================


def generate_references(
    dataset_dir: Path, muconeup_config: Path, log_file: Path
) -> dict[str, Path]:
    """Generate diploid references using dual simulation mode."""
    logger = logging.getLogger(__name__)
    logger.info("Generating diploid references (dual simulation mode)")

    refs_dir = dataset_dir / "references"
    ensure_dir(refs_dir)

    # Single command generates both variants (asymmetric diploid with two --fixed-lengths)
    cmd = [
        "muconeup",
        "--verbose",
        "--config",
        str(muconeup_config),
        "simulate",
        "--out-base",
        str(refs_dir / SAMPLE_BASE_NAME),
        "--fixed-lengths",
        str(VNTR_H1_LENGTH),
        "--fixed-lengths",
        str(VNTR_H2_LENGTH),
        "--output-structure",
        "--mutation-name",
        MUTATION_NAME,
        "--mutation-targets",
        MUTATION_TARGETS,
        "--seed",
        str(SEED_BASE),
    ]

    run_command(cmd, "Generate diploid references", log_file=log_file)

    # Organize files by variant
    references = {}
    variant_mapping = {"normal": "normal", "mut": "dupC"}

    for file_variant, dir_variant in variant_mapping.items():
        variant_dir = refs_dir / dir_variant
        ensure_dir(variant_dir)

        # Move files to variant subdirectory
        pattern = f"{SAMPLE_BASE_NAME}.001.{file_variant}.*"
        for file in refs_dir.glob(pattern):
            if file.is_file():
                dest = variant_dir / file.name
                file.rename(dest)
                logger.debug(f"  Moved: {file.name} -> {dir_variant}/")

                # Track FASTA file
                if file.suffix == ".fa":
                    references[dir_variant] = dest

    if len(references) != 2:
        raise RuntimeError(
            f"Expected 2 reference FASTAs, found {len(references)}: {list(references.keys())}"
        )

    logger.info(f"✓ Generated references: {list(references.keys())}")

    return references


def run_orf_analysis(fasta_path: Path, variant: str, muconeup_config: Path, log_file: Path):
    """Run ORF/toxicity analysis on a FASTA file."""
    logger = logging.getLogger(__name__)
    logger.info(f"  Running ORF analysis for {variant}")

    output_dir = fasta_path.parent
    base_name = fasta_path.stem.replace(".simulated", "")
    out_base = output_dir / f"{base_name}_toxic"

    cmd = [
        "muconeup",
        "--verbose",
        "--config",
        str(muconeup_config),
        "analyze",
        "orfs",
        str(fasta_path),
        "--out-base",
        str(out_base),
        "--orf-min-aa",
        "100",
        "--orf-aa-prefix",
        "MTSSV",
    ]

    run_command(cmd, f"ORF analysis for {variant}", log_file=log_file)

    # Log generated files
    for suffix in [".orf_stats.json", ".orfs.fa", ".orfs.faa"]:
        file_path = Path(str(out_base) + suffix)
        if file_path.exists():
            logger.debug(f"    Generated: {file_path.name}")


# ============================================================================
# Read Simulation
# ============================================================================


def simulate_platform_reads(
    platform_name: str,
    platform_config: dict,
    references: dict[str, Path],
    dataset_dir: Path,
    muconeup_config: Path,
    threads: int,
    log_file: Path,
):
    """Simulate reads for a platform across all variants."""
    logger = logging.getLogger(__name__)
    logger.info(
        f"Simulating {platform_name.upper()} reads (coverage={platform_config['coverage']}x)"
    )

    platform_dir = dataset_dir / platform_config["output_subdir"]
    ensure_dir(platform_dir)

    base_seed = SEED_BASE + platform_config["seed_offset"]

    for variant, ref_fa in references.items():
        variant_dir = platform_dir / variant
        ensure_dir(variant_dir)

        # Use different seeds for normal vs mutated (critical for avoiding BWA conflicts)
        seed = base_seed if variant == "normal" else base_seed + 10000
        logger.info(f"  Variant: {variant} (seed={seed})")

        cmd = [
            "muconeup",
            "--config",
            str(muconeup_config),
            "reads",
            platform_name,
            str(ref_fa),
            "--seed",
            str(seed),
        ]

        # Platform-specific parameters
        if platform_name == "illumina":
            cmd.extend(["--threads", str(threads)])
        elif platform_name in ["ont", "pacbio"]:
            cmd.extend(["--coverage", str(platform_config["coverage"])])

        run_command(cmd, f"Simulate {platform_name} reads for {variant}", log_file=log_file)

        # Collect and organize files
        collect_and_organize_files(ref_fa, variant_dir, platform_name)

        # Run SNaPshot validation for dupC variant
        if variant == "dupC":
            run_snapshot_validation(ref_fa, variant_dir, muconeup_config, log_file)


def collect_and_organize_files(ref_fa: Path, variant_dir: Path, platform: str):
    """Collect generated read files and organize into directories."""
    logger = logging.getLogger(__name__)

    # Extract base name
    base = ref_fa.stem.replace(".simulated", "")
    parent_dir = ref_fa.parent

    # Platform-specific file patterns and rename mappings
    if platform == "illumina":
        patterns = [
            f"{base}.simulated_R1.fastq.gz",
            f"{base}.simulated_R2.fastq.gz",
            f"{base}.simulated.bam",
            f"{base}.simulated.bam.bai",
            f"{base}.simulated_vntr_biased.bam",
            f"{base}.simulated_vntr_biased.bam.bai",
            f"{base}.simulated_vntr_biased_R1.fastq.gz",
            f"{base}.simulated_vntr_biased_R2.fastq.gz",
        ]
        rename_map = {
            f"{base}.simulated_R1.fastq.gz": "reads_R1.fastq.gz",
            f"{base}.simulated_R2.fastq.gz": "reads_R2.fastq.gz",
            f"{base}.simulated.bam": "aligned.bam",
            f"{base}.simulated.bam.bai": "aligned.bam.bai",
            f"{base}.simulated_vntr_biased.bam": "vntr_biased.bam",
            f"{base}.simulated_vntr_biased.bam.bai": "vntr_biased.bam.bai",
            f"{base}.simulated_vntr_biased_R1.fastq.gz": "vntr_biased_R1.fastq.gz",
            f"{base}.simulated_vntr_biased_R2.fastq.gz": "vntr_biased_R2.fastq.gz",
        }
    elif platform == "ont":
        # ONT files have different naming: *_ont.bam, *_ont_merged.fastq (NOT gzipped!)
        patterns = [
            f"{base}.simulated_ont_merged.fastq",
            f"{base}.simulated_ont.bam",
            f"{base}.simulated_ont.bam.bai",
            f"{base}.simulated_ont_metadata.tsv",
        ]
        rename_map = {
            f"{base}.simulated_ont_merged.fastq": "reads.fastq",
            f"{base}.simulated_ont.bam": "aligned.bam",
            f"{base}.simulated_ont.bam.bai": "aligned.bam.bai",
            f"{base}.simulated_ont_metadata.tsv": "metadata.tsv",
        }
    elif platform == "pacbio":
        # PacBio files: *_hifi_*.*, *_aligned.bam, *_metadata.tsv, *_vntr_efficiency_stats.json
        patterns = [
            f"{base}.simulated_hifi_0001.bam.pbi",
            f"{base}.simulated_hifi_0001.ccs_report.txt",
            f"{base}.simulated_hifi_0001.zmw_metrics.json.gz",
            f"{base}.simulated_hifi_0002.bam.pbi",
            f"{base}.simulated_hifi_0002.ccs_report.txt",
            f"{base}.simulated_hifi_0002.zmw_metrics.json.gz",
            f"{base}.simulated_aligned.bam",
            f"{base}.simulated_aligned.bam.bai",
            f"{base}.simulated_metadata.tsv",
            f"{base}.simulated_vntr_efficiency_stats.json",
        ]
        rename_map = {
            f"{base}.simulated_hifi_0001.bam.pbi": "hifi_0001.bam.pbi",
            f"{base}.simulated_hifi_0001.ccs_report.txt": "hifi_0001.ccs_report.txt",
            f"{base}.simulated_hifi_0001.zmw_metrics.json.gz": "hifi_0001.zmw_metrics.json.gz",
            f"{base}.simulated_hifi_0002.bam.pbi": "hifi_0002.bam.pbi",
            f"{base}.simulated_hifi_0002.ccs_report.txt": "hifi_0002.ccs_report.txt",
            f"{base}.simulated_hifi_0002.zmw_metrics.json.gz": "hifi_0002.zmw_metrics.json.gz",
            f"{base}.simulated_aligned.bam": "aligned.bam",
            f"{base}.simulated_aligned.bam.bai": "aligned.bam.bai",
            f"{base}.simulated_metadata.tsv": "metadata.tsv",
            f"{base}.simulated_vntr_efficiency_stats.json": "vntr_efficiency_stats.json",
        }
    else:
        raise ValueError(f"Unknown platform: {platform}")

    # Find and move files
    file_count = 0
    for pattern in patterns:
        source = parent_dir / pattern
        if source.exists():
            dest_name = rename_map.get(pattern, pattern)
            dest = variant_dir / dest_name
            source.rename(dest)
            file_count += 1
            logger.debug(f"      Moved: {source.name} -> {dest_name}")
        else:
            logger.warning(f"      Expected file not found: {pattern}")

    logger.info(f"    ✓ Generated {file_count} files")


def run_snapshot_validation(
    fasta_path: Path, output_dir: Path, muconeup_config: Path, log_file: Path
):
    """Run SNaPshot validation on a mutated FASTA file."""
    logger = logging.getLogger(__name__)
    logger.info("  Running SNaPshot validation for dupC")

    base_name = fasta_path.stem.replace(".simulated", "")
    output_file = output_dir / f"{base_name}_snapshot_validation.json"

    cmd = [
        "muconeup",
        "--verbose",
        "--config",
        str(muconeup_config),
        "analyze",
        "snapshot-validate",
        str(fasta_path),
        "--mutation",
        "dupC",
        "--output",
        str(output_file),
    ]

    try:
        run_command(cmd, "SNaPshot validation for dupC", check=False, log_file=log_file)

        if output_file.exists() and output_file.stat().st_size > 0:
            logger.info(f"    ✓ Validation output created: {output_file.name}")
        else:
            logger.warning("    SNaPshot validation skipped (non-critical)")

    except Exception as e:
        logger.warning(f"    SNaPshot validation failed (non-critical): {e}")


# ============================================================================
# Metadata Generation
# ============================================================================


def generate_manifest(dataset_dir: Path, all_files: list[Path], version: str) -> Path:
    """Generate manifest.json with complete metadata."""
    logger = logging.getLogger(__name__)
    logger.info("Generating manifest.json")

    manifest: dict[str, Any] = {
        "version": version,
        "generated_with": "MucOneUp",
        "muconeup_version": version,
        "generation_date": datetime.utcnow().isoformat() + "Z",
        "sample": {
            "base_name": SAMPLE_BASE_NAME,
            "vntr_h1": VNTR_H1_LENGTH,
            "vntr_h2": VNTR_H2_LENGTH,
            "mutation_mode": MUTATION_NAME,
            "mutation_targets": MUTATION_TARGETS,
            "seed_base": SEED_BASE,
        },
        "platforms": {},
        "structure": {
            "references/": "Shared diploid references (normal + dupC)",
            "illumina/": "Illumina paired-end sequencing (50x)",
            "ont/": "Oxford Nanopore sequencing (30x)",
            "pacbio/": "PacBio HiFi sequencing (30x)",
        },
        "checksums": {},
        "file_count": len(all_files),
    }

    # Platform metadata
    for platform_name, platform_config in PLATFORMS.items():
        if platform_config["enabled"]:
            manifest["platforms"][platform_name] = {
                "enabled": platform_config["enabled"],
                "coverage": platform_config["coverage"],
                "seed": SEED_BASE + cast(int, platform_config["seed_offset"]),
            }

    # Calculate checksums for all files
    for file_path in all_files:
        rel_path = file_path.relative_to(dataset_dir)
        checksum = calculate_checksum(file_path)
        manifest["checksums"][str(rel_path)] = checksum

    manifest_path = dataset_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    logger.info(f"✓ Manifest created: {len(manifest['checksums'])} files")

    return manifest_path


def generate_checksums_file(dataset_dir: Path, manifest: dict) -> Path:
    """Generate checksums.txt in sha256sum format."""
    logger = logging.getLogger(__name__)
    logger.info("Generating checksums.txt")

    checksums_path = dataset_dir / "checksums.txt"

    with open(checksums_path, "w") as f:
        for rel_path, checksum in sorted(manifest["checksums"].items()):
            f.write(f"{checksum}  {rel_path}\n")

    logger.info(f"✓ Checksums file created: {len(manifest['checksums'])} entries")

    return checksums_path


def generate_readme(dataset_dir: Path, version: str) -> Path:
    """Generate README.md for the dataset."""
    logger = logging.getLogger(__name__)
    logger.info("Generating README.md")

    h1 = VNTR_H1_LENGTH
    h2 = VNTR_H2_LENGTH

    readme_content = f"""# MucOneUp Test Dataset {version}

## Overview

Comprehensive test dataset with **asymmetric VNTR ({h1}/{h2} repeats)** across three sequencing platforms.

**Generated with**: MucOneUp {version}
**Platforms**: Illumina (50x), ONT (30x), PacBio (30x)
**Sample**: Asymmetric diploid (h1={h1}, h2={h2} repeats)
**Variants**: Normal baseline + dupC mutation

## Structure

```
testdata_40-70_{version}/
├── references/          # Shared diploid references
│   ├── normal/          # Normal baseline
│   └── dupC/            # dupC mutation (h1:r20)
├── illumina/            # Illumina paired-end (50x)
│   ├── normal/
│   └── dupC/
├── ont/                 # Oxford Nanopore (30x)
│   ├── normal/
│   └── dupC/
└── pacbio/              # PacBio HiFi (30x)
    ├── normal/
    └── dupC/
```

## Quick Start

```bash
# Verify integrity
sha256sum -c checksums.txt

# Explore references
cat references/normal/testdata_40-70.001.normal.simulated.structure.tsv
cat references/dupC/testdata_40-70.001.dupC.simulated.structure.tsv

# View Illumina reads
zcat illumina/normal/reads_R1.fastq.gz | head -4

# View ONT reads (uncompressed)
cat ont/normal/reads.fastq | head -4

# View PacBio reads
zcat pacbio/normal/hifi_0001.fastq.gz | head -4
```

## Platform Comparison

| Platform | Coverage | Read Type | VNTR-Biased |
|----------|----------|-----------|-------------|
| Illumina | 50x | Paired-end (2x150bp) | ✅ Yes |
| ONT | 30x | Long reads (~3kb) | ✅ Yes |
| PacBio | 30x | HiFi reads (~15kb) | ✅ Yes |

## Reproducibility

```bash
# Generate references (once)
muconeup --config config.json simulate \\
  --out-base testdata_40-70 \\
  --fixed-lengths 40 \\
  --fixed-lengths 70 \\
  --mutation-name normal,dupC \\
  --mutation-targets 1,20 \\
  --seed 42000

# Simulate Illumina (seeds: 42000 normal, 52000 dupC)
muconeup --config config.json reads illumina \\
  testdata_40-70.001.normal.simulated.fa --seed 42000
muconeup --config config.json reads illumina \\
  testdata_40-70.001.mut.simulated.fa --seed 52000

# Simulate ONT (seeds: 42010 normal, 52010 dupC)
muconeup --config config.json reads ont \\
  testdata_40-70.001.normal.simulated.fa --seed 42010
muconeup --config config.json reads ont \\
  testdata_40-70.001.mut.simulated.fa --seed 52010

# Simulate PacBio (seeds: 42020 normal, 52020 dupC)
muconeup --config config.json reads pacbio \\
  testdata_40-70.001.normal.simulated.fa --seed 42020
muconeup --config config.json reads pacbio \\
  testdata_40-70.001.mut.simulated.fa --seed 52020
```

## Citation

```bibtex
@software{{muconeup_{version.replace(".", "_").replace("-", "_")},
  author = {{Popp, Bernt}},
  title = {{MucOneUp: MUC1 VNTR Simulation Toolkit}},
  version = {{{version}}},
  year = {{2025}},
  url = {{https://github.com/berntpopp/MucOneUp}}
}}
```

## License

GPL-3.0 (same as MucOneUp)

## Support

- Documentation: https://berntpopp.github.io/MucOneUp/
- Issues: https://github.com/berntpopp/MucOneUp/issues
"""

    readme_path = dataset_dir / "README.md"
    with open(readme_path, "w") as f:
        f.write(readme_content)

    logger.info("✓ README created")

    return readme_path


# ============================================================================
# Tarball Creation
# ============================================================================


def create_tarball(dataset_dir: Path, version: str) -> Path:
    """Create compressed tarball of dataset in data/testdata_releases/ folder."""
    logger = logging.getLogger(__name__)
    logger.info("Creating tarball")

    # Create releases folder with date and version in data/
    date_str = datetime.now().strftime("%Y-%m-%d")
    releases_dir = Path("data/testdata_releases") / f"{date_str}_{version}"
    ensure_dir(releases_dir)

    tarball_name = f"{dataset_dir.name}.tar.gz"
    tarball_path = releases_dir / tarball_name

    cmd = ["tar", "-czf", str(tarball_path), "-C", str(dataset_dir.parent), dataset_dir.name]

    run_command(cmd, "Create tarball", check=True)

    logger.info(f"✓ Tarball created: {tarball_path}")

    # Create/update _latest directory with constant-name tarball
    latest_dir = Path("data/testdata_releases/_latest")
    if latest_dir.exists():
        shutil.rmtree(latest_dir)
    ensure_dir(latest_dir)

    # Copy tarball with constant name (version-agnostic for persistent URLs)
    constant_name = f"{SAMPLE_BASE_NAME}.tar.gz"
    latest_tarball = latest_dir / constant_name
    shutil.copy2(tarball_path, latest_tarball)

    # Add version metadata for traceability
    version_info = {
        "version": version,
        "release_date": date_str,
        "source_tarball": str(tarball_path.relative_to(Path.cwd())),
        "generated_utc": datetime.utcnow().isoformat() + "Z",
    }
    version_file = latest_dir / "VERSION.json"
    with open(version_file, "w") as f:
        json.dump(version_info, f, indent=2)

    logger.info("✓ Created _latest directory with constant-name tarball")
    logger.info(f"  Download URL: data/testdata_releases/_latest/{constant_name}")
    logger.info(f"  Version: {version} (released {date_str})")

    return tarball_path


# ============================================================================
# Argument Parsing
# ============================================================================


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate MucOneUp test dataset with multi-platform support"
    )
    parser.add_argument(
        "--version", default=DEFAULT_VERSION, help=f"Version tag (default: {DEFAULT_VERSION})"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output directory (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"MucOneUp config file (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help=f"Number of threads (default: {DEFAULT_THREADS})",
    )
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


# ============================================================================
# Main Execution
# ============================================================================


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

    # Filter platforms if specified
    enabled_platforms = {}
    if "all" in args.platforms:
        enabled_platforms = PLATFORMS.copy()
    else:
        for name in args.platforms:
            if name in PLATFORMS:
                enabled_platforms[name] = PLATFORMS[name]

    # Create output directory
    dataset_dir = args.output / f"testdata_40-70_{args.version}"
    ensure_dir(dataset_dir)

    # Create unified log file for all muconeup commands
    muconeup_log = dataset_dir / "muconeup.log"

    # Capture tool versions
    logger.info("Capturing tool versions...")
    capture_tool_versions(dataset_dir)

    logger.info("=" * 70)
    logger.info("MucOneUp Test Dataset Generator")
    logger.info("=" * 70)
    logger.info(f"Version: {args.version}")
    logger.info(f"Output: {dataset_dir}")
    logger.info(f"Platforms: {list(enabled_platforms.keys())}")
    logger.info(f"Threads: {args.threads}")
    logger.info("=" * 70)

    try:
        # Step 1: Generate diploid references
        logger.info("")
        logger.info("STEP 1: Generating diploid references")
        logger.info("-" * 70)
        references = generate_references(dataset_dir, args.config, muconeup_log)

        # Run ORF analysis on generated references
        logger.info("Running ORF/toxicity analysis on references")
        for dir_variant, ref_fa in references.items():
            try:
                run_orf_analysis(ref_fa, dir_variant, args.config, muconeup_log)
            except Exception as e:
                logger.warning(f"  ORF analysis failed for {dir_variant}: {e}")
        logger.info("✓ ORF analysis completed")

        # Step 2: Simulate reads for each platform
        logger.info("")
        logger.info("STEP 2: Simulating reads for enabled platforms")
        logger.info("-" * 70)

        for platform_name, platform_config in enabled_platforms.items():
            simulate_platform_reads(
                platform_name,
                platform_config,
                references,
                dataset_dir,
                args.config,
                args.threads,
                muconeup_log,
            )

        # Step 3: Generate metadata
        logger.info("")
        logger.info("STEP 3: Generating metadata")
        logger.info("-" * 70)

        # Collect all files
        all_files = collect_all_files(
            dataset_dir, exclude_patterns=["*.log", "*.pyc", "__pycache__"]
        )

        manifest_path = generate_manifest(dataset_dir, all_files, args.version)

        # Load manifest for checksums file
        with open(manifest_path) as f:
            manifest = json.load(f)

        generate_checksums_file(dataset_dir, manifest)
        generate_readme(dataset_dir, args.version)

        # Copy config for reproducibility
        shutil.copy(args.config, dataset_dir / "config.json")
        logger.info("✓ Copied config.json for reproducibility")

        # Step 4: Create tarball (optional)
        tarball_path = None
        if not args.no_tarball:
            logger.info("")
            logger.info("STEP 4: Creating tarball")
            logger.info("-" * 70)
            tarball_path = create_tarball(dataset_dir, args.version)

        # Cleanup temporary files
        if not args.no_cleanup:
            logger.info("")
            logger.info("Cleaning up temporary files (use --no-cleanup to keep them)")
            if tarball_path:
                import shutil as sh

                sh.rmtree(dataset_dir)
                logger.info(f"✓ Removed temporary directory: {dataset_dir}")

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
        logger.info(f"🧬 Platforms: {len(enabled_platforms)}")
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

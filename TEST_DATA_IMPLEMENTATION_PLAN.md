# Test Dataset Generation and Hosting Implementation Plan

**Issue**: #45 - Automated test data generation for documentation
**Goal**: Generate reproducible Illumina test dataset with asymmetric VNTR (40/70 repeats) using dual simulation (normal + dupC)
**Version**: v0.25.0+

---

## Executive Summary

This plan implements a **hybrid approach** combining GitHub Actions automation, GitHub Releases hosting, and documentation integration to provide a versioned, reproducible test dataset for MucOneUp users.

**Key Decisions** (based on 2025 best practices):
- ✅ GitHub Releases for primary hosting (2GB/file limit, unlimited total)
- ✅ Smaller files mirrored to `docs/assets/test-data/` for quick access
- ✅ Dedicated workflow on version tags (not docs build - keeps it fast)
- ✅ Deterministic generation with fixed seed
- ✅ SHA256 checksums for integrity validation
- ✅ Target size: ~50-80MB compressed
- ✅ **Dual simulation mode** (`normal,dupC`) generates both versions in one run

---

## Phase 1: Dataset Specification

### 1.1 Dataset Parameters

**Illumina Test Dataset** (Primary Focus)

```yaml
name: "illumina_testdata_40-70_dupC"
version: "v0.25.0"
description: "Illumina paired-end reads for asymmetric VNTR (40/70 repeats) with normal and dupC variants"

sample:
  name: "testdata_40-70"
  vntr_lengths: "40,70"  # Haplotype 1: 40 repeats, Haplotype 2: 70 repeats
  mutation_mode: "normal,dupC"  # Dual simulation - generates BOTH
  mutation_target: "1,20"  # dupC mutation on haplotype 1, repeat 20
  seed: 42000

  description: |
    Single diploid sample with asymmetric VNTR lengths:
    - Haplotype 1: 40 repeats (shorter VNTR)
    - Haplotype 2: 70 repeats (longer VNTR)

    Dual simulation generates:
    - testdata_40-70.001.normal.* (baseline, no mutation)
    - testdata_40-70.001.dupC.* (dupC mutation at h1:r20)

read_simulation:
  platform: "illumina"
  coverage: 50
  read_length: 150
  fragment_size_mean: 350
  fragment_size_sd: 50
  threads: 4

vntr_capture_efficiency:
  enabled: true
  penalty_factor: 0.375
  output_fastq:
    enabled: true
    preserve_read_names: true

output_files:
  # Normal variant files
  - "testdata_40-70.001.normal.simulated.fa"                     # Normal diploid reference
  - "testdata_40-70.001.normal.simulated.structure.tsv"          # VNTR structure (normal)
  - "testdata_40-70.001.normal.simulation_stats.json"            # Statistics (normal)
  - "testdata_40-70.001.normal.simulated_R1.fastq.gz"            # Normal reads R1
  - "testdata_40-70.001.normal.simulated_R2.fastq.gz"            # Normal reads R2
  - "testdata_40-70.001.normal.simulated.bam"                    # Normal aligned BAM
  - "testdata_40-70.001.normal.simulated.bam.bai"                # Normal BAM index
  - "testdata_40-70.001.normal.simulated_vntr_biased.bam"        # Normal VNTR-biased BAM
  - "testdata_40-70.001.normal.simulated_vntr_biased.bam.bai"   # Normal VNTR-biased index
  - "testdata_40-70.001.normal.simulated_vntr_biased_R1.fastq.gz"  # Normal VNTR-biased R1 ✨
  - "testdata_40-70.001.normal.simulated_vntr_biased_R2.fastq.gz"  # Normal VNTR-biased R2 ✨

  # dupC variant files (same structure)
  - "testdata_40-70.001.dupC.simulated.fa"
  - "testdata_40-70.001.dupC.simulated.structure.tsv"
  - "testdata_40-70.001.dupC.simulation_stats.json"
  - "testdata_40-70.001.dupC.simulated_R1.fastq.gz"
  - "testdata_40-70.001.dupC.simulated_R2.fastq.gz"
  - "testdata_40-70.001.dupC.simulated.bam"
  - "testdata_40-70.001.dupC.simulated.bam.bai"
  - "testdata_40-70.001.dupC.simulated_vntr_biased.bam"
  - "testdata_40-70.001.dupC.simulated_vntr_biased.bam.bai"
  - "testdata_40-70.001.dupC.simulated_vntr_biased_R1.fastq.gz"  # ✨
  - "testdata_40-70.001.dupC.simulated_vntr_biased_R2.fastq.gz"  # ✨

estimated_size: "60-80 MB compressed"
```

### 1.2 Why These Parameters?

**40/70 repeats (asymmetric)**: Realistic heterogeneous diploid scenario, one short + one long VNTR
**Dual simulation mode**: Efficient - generates normal AND dupC in one run with same seed
**dupC mutation @ h1:r20**: Most common MUC1-VNTR mutation, clinically relevant, on shorter haplotype
**Fixed seed 42000**: Reproducibility - same dataset every generation
**50x coverage**: Realistic for WGS, sufficient for variant calling
**VNTR bias enabled**: Demonstrates realistic coverage patterns (v0.25.0 feature)

**Why dual simulation?**
- ✅ More efficient (one command vs two)
- ✅ Perfect matched pair (same seed, same background)
- ✅ Demonstrates realistic comparison workflow
- ✅ Halves generation time (~5 minutes total instead of 10)

---

## Phase 2: Generation Script

### 2.1 Script Location

```
scripts/generate_test_data.py
```

### 2.2 Script Architecture

```python
#!/usr/bin/env python3
"""
Generate reproducible test dataset for MucOneUp documentation.

Usage:
    python scripts/generate_test_data.py --version v0.25.0 --output test_data/

This script:
1. Generates ONE diploid sample with asymmetric VNTR (40/70 repeats)
2. Uses dual simulation mode to create normal AND dupC variants
3. Runs Illumina read simulation for both variants
4. Creates metadata and checksums
5. Packages as versioned tarball
6. Validates output integrity
"""

import argparse
import hashlib
import json
import logging
import shutil
import subprocess
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List

# Sample configuration
SAMPLE_CONFIG = {
    "base_name": "testdata_40-70",
    "vntr_lengths": "40,70",  # Asymmetric: h1=40, h2=70
    "mutation_name": "normal,dupC",  # Dual simulation
    "mutation_targets": "1,20",  # dupC on haplotype 1, repeat 20
    "seed": 42000,
    "description": (
        "Diploid sample with asymmetric VNTR lengths (40/70 repeats). "
        "Dual simulation generates normal baseline and dupC mutant."
    ),
}


def run_command(cmd: List[str], description: str) -> None:
    """Run a command with logging."""
    logging.info(f"{description}...")
    logging.debug(f"Command: {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True
    )

    if result.stdout:
        logging.debug(f"stdout: {result.stdout}")
    if result.stderr:
        logging.debug(f"stderr: {result.stderr}")

    logging.info(f"✓ {description} completed")


def generate_haplotypes(config_path: str, base_name: str, output_dir: Path,
                       seed: int) -> None:
    """Generate diploid haplotypes using dual simulation mode."""

    cmd = [
        "muconeup",
        "--config", config_path,
        "simulate",
        "--out-base", str(output_dir / base_name),
        "--fixed-lengths", SAMPLE_CONFIG["vntr_lengths"],
        "--mutation-name", SAMPLE_CONFIG["mutation_name"],
        "--mutation-targets", SAMPLE_CONFIG["mutation_targets"],
        "--seed", str(seed),
    ]

    run_command(
        cmd,
        f"Generating haplotypes for {base_name} (dual simulation: normal + dupC)"
    )


def simulate_reads(config_path: str, fa_file: Path, seed: int,
                  threads: int = 4) -> None:
    """Simulate Illumina reads for a FASTA file."""

    variant_name = fa_file.stem.split('.')[-2]  # Extract 'normal' or 'dupC'

    cmd = [
        "muconeup",
        "--config", config_path,
        "reads", "illumina",
        str(fa_file),
        "--seed", str(seed),
        "--threads", str(threads),
    ]

    run_command(
        cmd,
        f"Simulating Illumina reads for {variant_name} variant"
    )


def calculate_checksum(file_path: Path) -> str:
    """Calculate SHA256 checksum for a file."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def collect_files_by_variant(output_dir: Path, base_name: str) -> Dict[str, List[Path]]:
    """Collect generated files organized by variant (normal/dupC)."""

    files_by_variant = {
        "normal": [],
        "dupC": []
    }

    # Pattern: testdata_40-70.001.{normal|dupC}.*
    for variant in ["normal", "dupC"]:
        pattern = f"{base_name}.001.{variant}.*"
        variant_files = sorted(output_dir.glob(pattern))
        files_by_variant[variant] = variant_files

        logging.info(f"Found {len(variant_files)} files for {variant} variant")

    return files_by_variant


def create_manifest(output_dir: Path, files_by_variant: Dict[str, List[Path]],
                   version: str) -> Path:
    """Create manifest.json with metadata and checksums."""

    manifest = {
        "version": version,
        "generated_with": "MucOneUp",
        "muconeup_version": version,
        "generation_date": datetime.utcnow().isoformat() + "Z",
        "description": SAMPLE_CONFIG["description"],
        "sample": {
            "base_name": SAMPLE_CONFIG["base_name"],
            "vntr_lengths": SAMPLE_CONFIG["vntr_lengths"],
            "vntr_h1": 40,
            "vntr_h2": 70,
            "mutation_mode": SAMPLE_CONFIG["mutation_name"],
            "mutation_target": SAMPLE_CONFIG["mutation_targets"],
            "seed": SAMPLE_CONFIG["seed"],
        },
        "variants": {},
        "checksums": {}
    }

    # Organize by variant
    for variant, files in files_by_variant.items():
        manifest["variants"][variant] = {
            "description": f"{'Normal baseline' if variant == 'normal' else 'dupC mutation at haplotype 1, repeat 20'}",
            "file_count": len(files),
            "files": [f.name for f in files]
        }

        # Calculate checksums
        for file in files:
            checksum = calculate_checksum(file)
            manifest["checksums"][file.name] = checksum
            logging.debug(f"  {file.name}: {checksum[:16]}...")

    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    logging.info(f"Created manifest with {len(manifest['checksums'])} file checksums")

    return manifest_path


def create_readme(output_dir: Path, version: str) -> Path:
    """Create README.md for the test dataset."""

    readme_content = f'''# MucOneUp Test Dataset {version}

## Overview

Reproducible Illumina test dataset with asymmetric VNTR (40/70 repeats) for MucOneUp validation and tutorials.

**Generated with**: MucOneUp {version}
**Generation date**: {datetime.utcnow().strftime("%Y-%m-%d")}
**Platform**: Illumina (paired-end, 2×150bp)
**Coverage**: 50×
**VNTR Capture Efficiency**: Enabled (penalty factor 0.375)

## Sample Description

**Base name**: `testdata_40-70`

This dataset contains ONE diploid sample with **asymmetric VNTR lengths**:
- **Haplotype 1**: 40 VNTR repeats (shorter)
- **Haplotype 2**: 70 VNTR repeats (longer)

This asymmetry represents a realistic heterogeneous diploid scenario.

## Variants

Using **dual simulation mode**, two variants were generated from the same diploid template:

| Variant | Description | Files |
|---------|-------------|-------|
| **normal** | Baseline (no mutation) | `*.001.normal.*` |
| **dupC** | dupC mutation at h1:r20 | `*.001.dupC.*` |

Both variants share:
- Same VNTR structure (40/70 asymmetry)
- Same seed (42000) for reproducibility
- Same sequencing parameters

The dupC mutation was introduced at **haplotype 1, repeat 20** (on the shorter VNTR).

## Files

Each variant includes 11 files:

```
testdata_40-70.001.[variant]/
├── .simulated.fa                          # Diploid reference FASTA
├── .simulated.structure.tsv               # VNTR structure
├── .simulation_stats.json                 # Simulation statistics
├── .simulated_R1.fastq.gz                 # Paired-end reads (R1)
├── .simulated_R2.fastq.gz                 # Paired-end reads (R2)
├── .simulated.bam                         # Aligned reads (BAM)
├── .simulated.bam.bai                     # BAM index
├── .simulated_vntr_biased.bam             # VNTR-corrected alignment
├── .simulated_vntr_biased.bam.bai        # VNTR-corrected BAM index
├── .simulated_vntr_biased_R1.fastq.gz    # VNTR-corrected FASTQ (R1) ✨
└── .simulated_vntr_biased_R2.fastq.gz    # VNTR-corrected FASTQ (R2) ✨
```

**Total**: 22 files (11 per variant)

!!! tip "VNTR-Biased Outputs"
    The `*_vntr_biased.*` files include realistic VNTR capture efficiency corrections,
    matching the ~3.5× VNTR:flanking coverage ratio observed in real Twist v2 data.

    **Use these files for realistic variant calling benchmarks!**

## Quick Start

```bash
# Extract dataset
tar -xzf illumina_testdata_40-70_dupC_{version}.tar.gz
cd illumina_testdata_{version}/

# Verify checksums
sha256sum -c checksums.txt

# Inspect VNTR structure (normal variant)
cat testdata_40-70.001.normal.simulated.structure.tsv
# Haplotype 1: 40 repeats
# Haplotype 2: 70 repeats

# Inspect VNTR structure (dupC variant)
cat testdata_40-70.001.dupC.simulated.structure.tsv
# Haplotype 1: 40 repeats with dupC mutation at repeat 20
# Haplotype 2: 70 repeats (unchanged)

# View reads (VNTR-biased)
samtools view testdata_40-70.001.normal.simulated_vntr_biased.bam | head
zcat testdata_40-70.001.dupC.simulated_vntr_biased_R1.fastq.gz | head

# Compare normal vs dupC coverage
samtools depth testdata_40-70.001.normal.simulated_vntr_biased.bam > normal_depth.txt
samtools depth testdata_40-70.001.dupC.simulated_vntr_biased.bam > dupC_depth.txt
```

## Use Cases

### 1. Variant Caller Validation

Compare variant calling results between normal and dupC:

```bash
# Call variants on normal sample
bcftools mpileup -f testdata_40-70.001.normal.simulated.fa \\
  testdata_40-70.001.normal.simulated_vntr_biased.bam | \\
  bcftools call -mv -Oz -o normal_variants.vcf.gz

# Call variants on dupC sample
bcftools mpileup -f testdata_40-70.001.dupC.simulated.fa \\
  testdata_40-70.001.dupC.simulated_vntr_biased.bam | \\
  bcftools call -mv -Oz -o dupC_variants.vcf.gz

# Expect to detect dupC mutation at h1:r20
bcftools view dupC_variants.vcf.gz | grep -i "dup"
```

### 2. Coverage Analysis

Analyze VNTR capture efficiency effect:

```bash
# Original coverage (no VNTR bias)
samtools depth testdata_40-70.001.normal.simulated.bam | \\
  awk '{{sum+=$3}} END {{print "Mean:", sum/NR}}'

# VNTR-biased coverage (realistic)
samtools depth testdata_40-70.001.normal.simulated_vntr_biased.bam | \\
  awk '{{sum+=$3}} END {{print "Mean:", sum/NR}}'
```

### 3. Mutation Detection

Identify dupC mutation location:

```bash
# Compare structure files
diff testdata_40-70.001.normal.simulated.structure.tsv \\
     testdata_40-70.001.dupC.simulated.structure.tsv

# Expected: Repeat 20 on haplotype 1 shows dupC mutation marker
```

### 4. Asymmetric VNTR Analysis

Study the effect of asymmetric VNTR lengths:

```bash
# Analyze haplotype-specific coverage
# H1 (40 repeats) vs H2 (70 repeats)

# Extract reads mapping to each haplotype
samtools view -h testdata_40-70.001.normal.simulated_vntr_biased.bam \\
  "haplotype1" > h1_reads.sam
samtools view -h testdata_40-70.001.normal.simulated_vntr_biased.bam \\
  "haplotype2" > h2_reads.sam
```

## Reproducibility

To regenerate this exact dataset:

```bash
# Single command (dual simulation)
muconeup --config config.json simulate \\
  --out-base testdata_40-70 \\
  --fixed-lengths 40,70 \\
  --mutation-name normal,dupC \\
  --mutation-targets 1,20 \\
  --seed 42000

# This generates both .001.normal.* AND .001.dupC.* files

# Then simulate reads for each variant
muconeup --config config.json reads illumina \\
  testdata_40-70.001.normal.simulated.fa \\
  --seed 42000

muconeup --config config.json reads illumina \\
  testdata_40-70.001.dupC.simulated.fa \\
  --seed 42000
```

## File Naming Convention

Files follow the pattern: `{base}.{series}.{variant}.{type}`

- **base**: `testdata_40-70` (indicates 40/70 asymmetric VNTR)
- **series**: `001` (first in series)
- **variant**: `normal` or `dupC`
- **type**: File type (`.simulated.fa`, `_R1.fastq.gz`, etc.)

## Statistics

Each variant includes `simulation_stats.json` with:
- Total sequence length
- VNTR region coordinates
- Repeat counts per haplotype
- Mutation information (dupC only)

```bash
# View normal stats
jq . testdata_40-70.001.normal.simulation_stats.json

# View dupC stats (includes mutation metadata)
jq . testdata_40-70.001.dupC.simulation_stats.json
```

## Citation

If you use this dataset, please cite MucOneUp:

```bibtex
@software{{muconeup_{version.replace('.', '_').replace('-', '_')},
  author = {{Popp, Bernt}},
  title = {{MucOneUp: MUC1 VNTR Simulation Toolkit}},
  version = {{{version}}},
  year = {{2025}},
  url = {{https://github.com/berntpopp/MucOneUp}}
}}
```

## License

This test dataset is provided under the same license as MucOneUp (GPL-3.0).

## Support

- **Documentation**: https://berntpopp.github.io/MucOneUp/
- **Issues**: https://github.com/berntpopp/MucOneUp/issues
- **Discussions**: https://github.com/berntpopp/MucOneUp/discussions

---

**Generated by**: `scripts/generate_test_data.py`
**Seed**: 42000 (for reproducibility)
**Dual simulation**: normal,dupC modes
'''

    readme_path = output_dir / "README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)

    logging.info("Created README.md")

    return readme_path


def create_checksums_file(output_dir: Path, checksums: Dict[str, str]) -> Path:
    """Create checksums.txt in sha256sum format."""

    checksums_path = output_dir / "checksums.txt"

    with open(checksums_path, 'w') as f:
        for filename, checksum in sorted(checksums.items()):
            f.write(f"{checksum}  {filename}\n")

    logging.info(f"Created checksums.txt with {len(checksums)} entries")

    return checksums_path


def create_tarball(output_dir: Path, version: str) -> Path:
    """Create compressed tarball of test dataset."""

    tarball_name = f"illumina_testdata_40-70_dupC_{version}.tar.gz"
    tarball_path = output_dir.parent / tarball_name

    logging.info(f"Creating tarball: {tarball_name}")

    with tarfile.open(tarball_path, "w:gz", compresslevel=9) as tar:
        tar.add(output_dir, arcname=output_dir.name)

    # Calculate tarball checksum
    tarball_checksum = calculate_checksum(tarball_path)

    checksum_file = tarball_path.with_suffix('.tar.gz.sha256')
    with open(checksum_file, 'w') as f:
        f.write(f"{tarball_checksum}  {tarball_name}\n")

    size_mb = tarball_path.stat().st_size / 1024 / 1024
    logging.info(f"✓ Tarball created: {size_mb:.2f} MB")
    logging.info(f"✓ Tarball SHA256: {tarball_checksum}")

    return tarball_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate MucOneUp test dataset with asymmetric VNTR (40/70 repeats)"
    )
    parser.add_argument(
        "--version",
        default="v0.25.0",
        help="Version tag (default: v0.25.0)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("test_data"),
        help="Output directory (default: test_data/)"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.json"),
        help="MucOneUp config file (default: config.json)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose debug logging"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Create output directory
    dataset_dir = args.output / f"illumina_testdata_{args.version}"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"=" * 60)
    logging.info(f"MucOneUp Test Dataset Generator")
    logging.info(f"=" * 60)
    logging.info(f"Version: {args.version}")
    logging.info(f"Output: {dataset_dir}")
    logging.info(f"Config: {args.config}")
    logging.info(f"Threads: {args.threads}")
    logging.info(f"=" * 60)

    # Step 1: Generate haplotypes (dual simulation: normal + dupC)
    logging.info("STEP 1: Generating haplotypes (dual simulation)")
    generate_haplotypes(
        config_path=str(args.config),
        base_name=SAMPLE_CONFIG["base_name"],
        output_dir=dataset_dir,
        seed=SAMPLE_CONFIG["seed"]
    )

    # Step 2: Simulate reads for each variant
    logging.info("STEP 2: Simulating Illumina reads")

    base_name = SAMPLE_CONFIG["base_name"]

    # Find generated FASTA files
    normal_fa = dataset_dir / f"{base_name}.001.normal.simulated.fa"
    dupC_fa = dataset_dir / f"{base_name}.001.dupC.simulated.fa"

    if not normal_fa.exists():
        raise FileNotFoundError(f"Normal FASTA not found: {normal_fa}")
    if not dupC_fa.exists():
        raise FileNotFoundError(f"dupC FASTA not found: {dupC_fa}")

    # Simulate reads for normal variant
    simulate_reads(
        config_path=str(args.config),
        fa_file=normal_fa,
        seed=SAMPLE_CONFIG["seed"],
        threads=args.threads
    )

    # Simulate reads for dupC variant
    simulate_reads(
        config_path=str(args.config),
        fa_file=dupC_fa,
        seed=SAMPLE_CONFIG["seed"],
        threads=args.threads
    )

    # Step 3: Collect and organize files
    logging.info("STEP 3: Collecting generated files")
    files_by_variant = collect_files_by_variant(dataset_dir, base_name)

    total_files = sum(len(files) for files in files_by_variant.values())
    logging.info(f"Total files generated: {total_files}")

    # Step 4: Create metadata
    logging.info("STEP 4: Creating manifest and checksums")
    manifest_path = create_manifest(dataset_dir, files_by_variant, args.version)

    with open(manifest_path) as f:
        manifest = json.load(f)

    create_checksums_file(dataset_dir, manifest["checksums"])
    create_readme(dataset_dir, args.version)

    # Copy config for reproducibility
    shutil.copy(args.config, dataset_dir / "config.json")
    logging.info("Copied config.json for reproducibility")

    # Step 5: Create tarball
    logging.info("STEP 5: Creating tarball")
    tarball_path = create_tarball(dataset_dir, args.version)

    # Final summary
    logging.info("=" * 60)
    logging.info("✅ Test dataset generation complete!")
    logging.info("=" * 60)
    logging.info(f"📦 Tarball: {tarball_path.name}")
    logging.info(f"   Size: {tarball_path.stat().st_size / 1024 / 1024:.2f} MB")
    logging.info(f"📋 Manifest: {manifest_path.name}")
    logging.info(f"   Files: {total_files}")
    logging.info(f"   Variants: {len(files_by_variant)}")
    logging.info(f"🔐 Checksums: checksums.txt")
    logging.info(f"   Entries: {len(manifest['checksums'])}")
    logging.info("=" * 60)
    logging.info(f"Next: Upload {tarball_path.name} to GitHub Release {args.version}")


if __name__ == "__main__":
    main()
```

---

## Phase 3: GitHub Actions Workflow

### 3.1 Workflow File

```yaml
# .github/workflows/generate-test-data.yml
name: Generate Test Data

on:
  push:
    tags:
      - 'v*.*.*'  # Trigger on version tags
  workflow_dispatch:  # Manual trigger
    inputs:
      version:
        description: 'Version tag (e.g., v0.25.0)'
        required: true
        default: 'v0.25.0'

permissions:
  contents: write  # Required for release uploads

jobs:
  generate-illumina-testdata:
    name: Generate Illumina Test Dataset
    runs-on: ubuntu-latest
    timeout-minutes: 30

    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.10'

      - name: Install uv
        uses: astral-sh/setup-uv@v5
        with:
          version: "latest"

      - name: Install MucOneUp
        run: |
          uv pip install --system .

      - name: Install external dependencies
        run: |
          # Illumina pipeline dependencies
          sudo apt-get update
          sudo apt-get install -y samtools bwa

      - name: Verify installation
        run: |
          muconeup --version
          samtools --version
          bwa version

      - name: Generate test dataset
        run: |
          python scripts/generate_test_data.py \
            --version ${{ github.ref_name || inputs.version }} \
            --config config.json \
            --output test_data/ \
            --threads 4 \
            --verbose

      - name: Verify dataset integrity
        run: |
          cd test_data/illumina_testdata_${{ github.ref_name || inputs.version }}/
          sha256sum -c checksums.txt

      - name: Show dataset summary
        run: |
          cd test_data/illumina_testdata_${{ github.ref_name || inputs.version }}/
          echo "=== Dataset Summary ==="
          cat manifest.json | jq '.sample, .variants'
          echo "=== File Counts ==="
          ls -1 | wc -l
          echo "=== Total Size ==="
          du -sh .

      - name: Upload artifacts (for inspection)
        uses: actions/upload-artifact@v4
        with:
          name: illumina-testdata-${{ github.ref_name || inputs.version }}
          path: |
            test_data/*.tar.gz
            test_data/*.tar.gz.sha256
          retention-days: 30

      - name: Create GitHub Release
        if: startsWith(github.ref, 'refs/tags/')
        uses: softprops/action-gh-release@v1
        with:
          files: |
            test_data/*.tar.gz
            test_data/*.tar.gz.sha256
          body: |
            ## Test Dataset for MucOneUp ${{ github.ref_name }}

            Reproducible Illumina test dataset with **asymmetric VNTR (40/70 repeats)** for validation and tutorials.

            ### Download

            **Illumina Dataset** (`illumina_testdata_40-70_dupC_${{ github.ref_name }}.tar.gz`) - ~60-80 MB

            - **ONE diploid sample** with asymmetric VNTR:
              - Haplotype 1: 40 repeats
              - Haplotype 2: 70 repeats
            - **TWO variants** (dual simulation):
              - `normal`: Baseline (no mutation)
              - `dupC`: dupC mutation at h1:r20
            - **50× coverage** with VNTR capture efficiency enabled
            - **Includes**: FASTQ, BAM, and VNTR-biased outputs

            ### Quick Start

            ```bash
            # Download and extract
            wget https://github.com/berntpopp/MucOneUp/releases/download/${{ github.ref_name }}/illumina_testdata_40-70_dupC_${{ github.ref_name }}.tar.gz
            tar -xzf illumina_testdata_40-70_dupC_${{ github.ref_name }}.tar.gz

            # Verify integrity
            cd illumina_testdata_${{ github.ref_name }}/
            sha256sum -c checksums.txt

            # View VNTR structure
            cat testdata_40-70.001.normal.simulated.structure.tsv
            cat testdata_40-70.001.dupC.simulated.structure.tsv
            ```

            ### Documentation

            See [Test Data Guide](https://berntpopp.github.io/MucOneUp/guides/test-data/) for detailed usage instructions.

            ### What's Included

            **22 files total** (11 per variant):
            - Diploid reference FASTA
            - VNTR structure TSV
            - Simulation statistics JSON
            - Paired-end FASTQ (R1/R2)
            - Aligned BAM + index
            - VNTR-biased BAM + index
            - VNTR-biased FASTQ (R1/R2) ✨ **NEW in v0.25.0**

            ### Reproducibility

            ```bash
            # Single command generates both variants
            muconeup --config config.json simulate \
              --out-base testdata_40-70 \
              --fixed-lengths 40,70 \
              --mutation-name normal,dupC \
              --mutation-targets 1,20 \
              --seed 42000
            ```
          draft: false
          prerelease: false
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}

      - name: Copy to docs assets (for quick access)
        if: startsWith(github.ref, 'refs/tags/')
        run: |
          mkdir -p docs/assets/test-data/
          cp test_data/*.tar.gz docs/assets/test-data/
          cp test_data/*.tar.gz.sha256 docs/assets/test-data/

          # Create version symlink (latest)
          cd docs/assets/test-data/
          ln -sf illumina_testdata_40-70_dupC_${{ github.ref_name }}.tar.gz \
                 illumina_testdata_40-70_dupC_latest.tar.gz

      - name: Commit docs assets
        if: startsWith(github.ref, 'refs/tags/')
        run: |
          git config --local user.email "github-actions[bot]@users.noreply.github.com"
          git config --local user.name "github-actions[bot]"
          git add docs/assets/test-data/
          git commit -m "docs: Add test data for ${{ github.ref_name }}" || echo "No changes to commit"
          git push origin HEAD:main
```

---

## Phase 4: Documentation Integration

### 4.1 New Guide Page

**File**: `docs/guides/test-data.md`

```markdown
# Test Data

MucOneUp provides a pre-generated test dataset with **asymmetric VNTR (40/70 repeats)** for quick validation, tutorials, and variant calling benchmarks.

## Quick Start

=== "Download from GitHub Releases"

    ```bash
    # Latest version
    wget https://github.com/berntpopp/MucOneUp/releases/latest/download/illumina_testdata_40-70_dupC_latest.tar.gz

    # Specific version (v0.25.0)
    wget https://github.com/berntpopp/MucOneUp/releases/download/v0.25.0/illumina_testdata_40-70_dupC_v0.25.0.tar.gz

    # Extract
    tar -xzf illumina_testdata_40-70_dupC_v0.25.0.tar.gz
    cd illumina_testdata_v0.25.0/
    ```

=== "Download from Docs (Mirror)"

    ```bash
    # Quick access
    wget https://berntpopp.github.io/MucOneUp/assets/test-data/illumina_testdata_40-70_dupC_latest.tar.gz

    # Extract
    tar -xzf illumina_testdata_40-70_dupC_latest.tar.gz
    ```

## Dataset Overview

### Illumina Test Dataset

**Version**: v0.25.0+
**Platform**: Illumina paired-end (2×150bp)
**Coverage**: 50× with VNTR capture efficiency correction
**Size**: ~60-80 MB compressed

### Sample Structure

**ONE diploid sample** with **asymmetric VNTR lengths**:

| Haplotype | VNTR Repeats | Description |
|-----------|--------------|-------------|
| Haplotype 1 | 40 | Shorter VNTR |
| Haplotype 2 | 70 | Longer VNTR |

This asymmetry represents a realistic heterogeneous diploid scenario.

### Variants

Using **dual simulation mode**, two variants were generated:

| Variant | Mutation | Seed | Description |
|---------|----------|------|-------------|
| `normal` | None | 42000 | Baseline (no mutation) |
| `dupC` | dupC @ h1:r20 | 42000 | dupC mutation on haplotype 1, repeat 20 |

Both variants:
- Share the same diploid structure (40/70 asymmetry)
- Use the same seed for reproducibility
- Were generated in a single command using `--mutation-name normal,dupC`

### Files Included

**22 files total** (11 per variant):

```
testdata_40-70.001.normal.*         # Normal variant (11 files)
testdata_40-70.001.dupC.*           # dupC variant (11 files)
```

Each variant contains:

```
├── .simulated.fa                          # Diploid reference FASTA
├── .simulated.structure.tsv               # VNTR structure
├── .simulation_stats.json                 # Simulation statistics
├── .simulated_R1.fastq.gz                 # Paired-end reads (R1)
├── .simulated_R2.fastq.gz                 # Paired-end reads (R2)
├── .simulated.bam                         # Aligned reads (BAM)
├── .simulated.bam.bai                     # BAM index
├── .simulated_vntr_biased.bam             # VNTR-corrected alignment
├── .simulated_vntr_biased.bam.bai        # VNTR-corrected BAM index
├── .simulated_vntr_biased_R1.fastq.gz    # VNTR-corrected FASTQ (R1) ✨
└── .simulated_vntr_biased_R2.fastq.gz    # VNTR-corrected FASTQ (R2) ✨
```

!!! tip "VNTR-Biased Outputs"
    The `*_vntr_biased.*` files include realistic VNTR capture efficiency corrections,
    matching the ~3.5× VNTR:flanking coverage ratio observed in real Twist v2 data.

    **Use these files for realistic variant calling benchmarks!**

## Verification

Always verify dataset integrity after download:

```bash
# Verify checksums
sha256sum -c checksums.txt

# Expected: All files report "OK"
```

## Usage Examples

### Example 1: Inspect VNTR Structure

```bash
# View normal VNTR structure (40/70 asymmetry)
cat testdata_40-70.001.normal.simulated.structure.tsv

# View dupC VNTR structure (mutation at h1:r20)
cat testdata_40-70.001.dupC.simulated.structure.tsv

# Compare structures
diff testdata_40-70.001.normal.simulated.structure.tsv \
     testdata_40-70.001.dupC.simulated.structure.tsv
```

### Example 2: Coverage Analysis

```bash
# Check coverage with samtools
samtools depth testdata_40-70.001.normal.simulated_vntr_biased.bam | head -20

# Mean coverage (normal variant)
samtools depth testdata_40-70.001.normal.simulated_vntr_biased.bam | \
  awk '{sum+=$3; count++} END {print "Mean coverage:", sum/count}'

# Compare VNTR-biased vs original coverage
samtools depth testdata_40-70.001.normal.simulated.bam | \
  awk '{sum+=$3} END {print "Original:", sum/NR}'
samtools depth testdata_40-70.001.normal.simulated_vntr_biased.bam | \
  awk '{sum+=$3} END {print "VNTR-biased:", sum/NR}'
```

### Example 3: Variant Calling (Normal vs dupC)

```bash
# Call variants on normal sample
bcftools mpileup -f testdata_40-70.001.normal.simulated.fa \
  testdata_40-70.001.normal.simulated_vntr_biased.bam | \
  bcftools call -mv -Oz -o normal_variants.vcf.gz

# Call variants on dupC sample
bcftools mpileup -f testdata_40-70.001.dupC.simulated.fa \
  testdata_40-70.001.dupC.simulated_vntr_biased.bam | \
  bcftools call -mv -Oz -o dupC_variants.vcf.gz

# Compare variant calls
bcftools isec normal_variants.vcf.gz dupC_variants.vcf.gz -p comparison/

# Inspect dupC-specific variants
bcftools view comparison/0001.vcf.gz | less
```

### Example 4: Read Quality Assessment

```bash
# FastQC on VNTR-biased reads (both variants)
fastqc testdata_40-70.001.normal.simulated_vntr_biased_R1.fastq.gz \
       testdata_40-70.001.normal.simulated_vntr_biased_R2.fastq.gz \
       testdata_40-70.001.dupC.simulated_vntr_biased_R1.fastq.gz \
       testdata_40-70.001.dupC.simulated_vntr_biased_R2.fastq.gz

# MultiQC summary
multiqc .
```

### Example 5: Mutation Detection

```bash
# Extract reads covering mutation site (h1:r20)
# (Requires knowing VNTR coordinates from structure.tsv)

# View simulation statistics
jq . testdata_40-70.001.dupC.simulation_stats.json

# Look for mutation information in stats
jq '.mutations' testdata_40-70.001.dupC.simulation_stats.json
```

## Reproducing the Dataset

To regenerate this exact dataset:

```bash
# SINGLE COMMAND generates BOTH normal AND dupC variants
muconeup --config config.json simulate \
  --out-base testdata_40-70 \
  --fixed-lengths 40,70 \
  --mutation-name normal,dupC \
  --mutation-targets 1,20 \
  --seed 42000

# This creates:
#   testdata_40-70.001.normal.simulated.fa
#   testdata_40-70.001.dupC.simulated.fa

# Then simulate reads for each variant
muconeup --config config.json reads illumina \
  testdata_40-70.001.normal.simulated.fa \
  --seed 42000

muconeup --config config.json reads illumina \
  testdata_40-70.001.dupC.simulated.fa \
  --seed 42000
```

!!! tip "Dual Simulation Mode"
    The `--mutation-name normal,dupC` flag uses **dual simulation mode**, which:

    - Generates both variants in one command (more efficient)
    - Ensures perfect matched pair (same seed, same background)
    - Halves generation time (~5 minutes instead of 10)

## Version History

| Version | Date | Changes |
|---------|------|---------|
| v0.25.0 | 2025-10-27 | Initial release with asymmetric VNTR (40/70) and dual simulation |

## Support

- **Documentation**: [https://berntpopp.github.io/MucOneUp/](https://berntpopp.github.io/MucOneUp/)
- **Issues**: [GitHub Issues](https://github.com/berntpopp/MucOneUp/issues)
- **Discussions**: [GitHub Discussions](https://github.com/berntpopp/MucOneUp/discussions)

## Citation

```bibtex
@software{muconeup_v0_25_0,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation Toolkit},
  version = {v0.25.0},
  year = {2025},
  url = {https://github.com/berntpopp/MucOneUp}
}
```
```

### 4.2 Update Quick Start

**File**: `docs/getting-started/quickstart.md`

Add at the beginning:

```markdown
## Option 1: Quick Start with Test Data ⚡ (Recommended)

!!! success "Fastest way to get started!"
    Skip generation and use pre-built test data to validate your installation in < 2 minutes.

```bash
# Download test dataset (~70 MB)
wget https://github.com/berntpopp/MucOneUp/releases/latest/download/illumina_testdata_40-70_dupC_latest.tar.gz

# Extract
tar -xzf illumina_testdata_40-70_dupC_latest.tar.gz
cd illumina_testdata_v0.25.0/

# Verify integrity
sha256sum -c checksums.txt

# Explore: ONE sample with 40/70 asymmetric VNTR, TWO variants (normal + dupC)
ls -lh testdata_40-70.001.normal.*
ls -lh testdata_40-70.001.dupC.*

# View VNTR structure
cat testdata_40-70.001.normal.simulated.structure.tsv
```

**What you get**:
- Asymmetric diploid (40/70 VNTR repeats)
- Normal and dupC variants
- FASTQ, BAM, and VNTR-biased outputs
- 50× coverage with realistic VNTR capture efficiency

**Next steps**: See [Test Data Guide](../guides/test-data.md) for usage examples.

---

## Option 2: Generate Your Own Data

If you want to generate custom data...
```

---

## Phase 5: Implementation Summary

### Key Changes from Original Plan

✅ **ONE sample, not four**: Single diploid with asymmetric 40/70 VNTR
✅ **Dual simulation**: `--mutation-name normal,dupC` generates both variants
✅ **More efficient**: ~5 minutes total instead of 10-15 minutes
✅ **Realistic scenario**: Asymmetric VNTR represents heterogeneous diploid
✅ **Perfect matched pair**: Same seed, same background, only mutation differs
✅ **22 files total**: 11 per variant (normal + dupC)

### Commands

**Single generation command**:
```bash
muconeup --config config.json simulate \
  --out-base testdata_40-70 \
  --fixed-lengths 40,70 \
  --mutation-name normal,dupC \
  --mutation-targets 1,20 \
  --seed 42000
```

**Then read simulation (2 commands, one per variant)**:
```bash
muconeup --config config.json reads illumina \
  testdata_40-70.001.normal.simulated.fa --seed 42000

muconeup --config config.json reads illumina \
  testdata_40-70.001.dupC.simulated.fa --seed 42000
```

### Timeline: ~6-8 Hours Total

Same as before, but generation is faster (~5 min instead of 15 min).

---

## References

Same as original plan.

---

**End of Revised Implementation Plan**

**Next Steps**:
1. Review revised plan
2. Create feature branch: `feature/test-data-generation`
3. Implement script with dual simulation mode
4. Test locally with: `python scripts/generate_test_data.py --version v0.25.0`
5. Create PR and deploy

**Questions?** Open discussion in Issue #45.

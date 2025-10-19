# Issue #28: Comprehensive Reference Assembly Management - Updated Assessment

**Author**: Expert Senior Developer & Bioinformatician Analysis
**Date**: 2025-10-19
**Status**: Design Complete - Ready for Implementation

---

## Executive Summary

This document provides a comprehensive analysis and solution design for implementing robust reference assembly management in MucOneUp, based on:
- Current industry best practices (GATK, Refgenie, NCBI/UCSC standards)
- Comparative analysis with VNtyper's implementation
- MucOneUp's existing architecture
- Extensibility requirements for future assemblies and organisms

**Key Recommendation**: Implement a lightweight, config-driven reference management system that follows the **Refgenie philosophy** while maintaining MucOneUp's simplicity and avoiding external dependencies.

---

## Table of Contents

1. [Current State Analysis](#1-current-state-analysis)
2. [Industry Best Practices](#2-industry-best-practices)
3. [VNtyper Comparison](#3-vntyper-comparison)
4. [Critical Requirements](#4-critical-requirements)
5. [Comprehensive Solution Design](#5-comprehensive-solution-design)
6. [Implementation Plan](#6-implementation-plan)
7. [Migration Strategy](#7-migration-strategy)
8. [Testing Strategy](#8-testing-strategy)
9. [Documentation Requirements](#9-documentation-requirements)

---

## 1. Current State Analysis

### 1.1 Existing Infrastructure

**Config Support** (`config.json`):
```json
{
  "constants": {
    "hg38": {
      "left": "...",
      "right": "...",
      "vntr_region": "chr1:155188487-155192239"
    },
    "hg19": {
      "left": "...",
      "right": "...",
      "vntr_region": "chr1:155160963-155162030"
    }
  },
  "reference_assembly": "hg38"  // Active assembly
}
```

**Validation Module** (`muc_one_up/bioinformatics/reference_validation.py`):
- ✅ `validate_reference_genome()` - Checks FASTA + indices
- ✅ `_validate_bwa_indices()` - BWA index validation
- ✅ `_validate_minimap2_indices()` - Minimap2 index validation
- ✅ `validate_bam_file()` - BAM/BAI validation
- ✅ `validate_bed_file()` - BED format validation

**Reference Utils** (`muc_one_up/read_simulator/utils/reference_utils.py`):
- ✅ `ReferenceInfo` - FASTA metadata (sequence count, IDs, lengths)
- ✅ `get_reference_info()` - Parse reference information
- ✅ `is_diploid_reference()` - Diploid detection
- ✅ `extract_haplotypes()` - Split diploid references
- ✅ `validate_reference_compatibility()` - Sequence count validation

**Helper Script** (`helpers/download_references.py`):
- ✅ Downloads flanking regions via UCSC DAS API
- ✅ Reverse complement handling (MUC1 on negative strand)
- ✅ NanoSim model setup
- ✅ Minimap2 index creation
- ⚠️ **Missing**: Explicit assembly-to-FASTA path mapping

### 1.2 Current Gaps

1. **No Centralized Reference Path Management**
   - Reference FASTA paths not linked to assemblies in config
   - Pipelines rely on external `human_reference` parameter
   - No validation that reference matches assembly

2. **Limited Assembly Support**
   - Hardcoded hg19/hg38 only
   - No support for other organisms (mouse, etc.)
   - No mechanism for custom assemblies

3. **Chromosome Naming Ambiguity**
   - No handling of chr1 (UCSC) vs 1 (NCBI) differences
   - No mitochondrial sequence handling (chrM vs chrMT vs MT)

4. **No Reference Provenance Tracking**
   - Unknown reference source/version
   - No checksums or validation
   - Difficult to reproduce analyses

---

## 2. Industry Best Practices

### 2.1 GATK Recommendations (2024)

**Reference Selection**:
- ✅ **Strongly recommend GRCh38/hg38** for new projects
- GRCh38 corrects thousands of SNP/indel artifacts present in GRCh37
- Includes synthetic centromeric sequences and updated non-nuclear genomic sequence

**Critical Principles**:
1. **Consistency is paramount** - Never mix reference assemblies mid-analysis
2. **Validate reference identity** - Use sequence dictionaries (.dict files)
3. **Consider analysis sets** - Account for ALT contigs, decoy sequences
4. **Chromosome naming matters** - Tools must handle chr1 vs 1 consistently

### 2.2 Refgenie Philosophy

**Key Concepts**:
1. **Asset-Based Management** - References are organized as versioned assets
2. **Pull vs Build** - Pre-built assets vs local construction
3. **Unique Asset IDs** - Computed from content (reproducibility)
4. **Parent-Child Tracking** - Record which assets derive from others
5. **Compatibility Levels** - Establish genome identity and compatibility

**Refgenie Asset Structure**:
```
reference_name/
  ├── asset_name/
  │   ├── fasta
  │   ├── fasta.fai
  │   ├── bwa_index/
  │   ├── minimap2_index/
  │   └── metadata.yaml
```

### 2.3 Chromosome Naming Conventions

| Source       | Autosome | X/Y  | Mitochondrial | Prefix  |
|------------- |----------|------|---------------|---------|
| **NCBI**     | 1-22     | X, Y | MT            | None    |
| **Ensembl**  | 1-22     | X, Y | MT            | None    |
| **UCSC**     | 1-22     | X, Y | M or MT       | chr     |
| **GRCh37**   | 1-22     | X, Y | MT            | None    |
| **hg19**     | 1-22     | X, Y | M (+ chrMT)   | chr     |
| **GRCh38**   | 1-22     | X, Y | MT            | None    |
| **hg38**     | 1-22     | X, Y | M             | chr     |

**MUC1 Locus** (chr1):
- UCSC: `chr1:155188487-155192239` (hg38), `chr1:155160963-155162030` (hg19)
- NCBI: `1:155188487-155192239` (GRCh38), `1:155160963-155162030` (GRCh37)

---

## 3. VNtyper Comparison

### 3.1 VNtyper's Approach

**Config Structure** (`vntyper/config.json`):
```json
{
  "reference_data": {
    "advntr_reference_vntr_hg19": "reference/vntr_data/.../hg19_genic_VNTRs.db",
    "advntr_reference_vntr_hg38": "reference/vntr_data/.../hg38_selected_VNTRs_Illumina.db",
    "bwa_reference_hg19": "reference/chr1.hg19.fa.gz",
    "bwa_reference_hg38": "reference/chr1.hg38.fa.gz"
  },
  "bam_processing": {
    "bam_region_hg19": "chr1:155158000-155163000",
    "bam_region_hg38": "chr1:155184000-155194000",
    "vntr_region_hg19": "chr1:155160500-155162000",
    "vntr_region_hg38": "chr1:155188000-155192500"
  }
}
```

**CLI Integration**:
```python
parser.add_argument(
    '--reference-assembly',
    type=str,
    choices=["hg19", "hg38"],
    default="hg19",
    help="Specify the reference assembly..."
)
```

**Runtime Selection**:
```python
# Pipeline determines paths based on assembly
if args.reference_assembly == "hg19":
    bwa_reference = config["reference_data"]["bwa_reference_hg19"]
    vntr_region = config["bam_processing"]["vntr_region_hg19"]
else:
    bwa_reference = config["reference_data"]["bwa_reference_hg38"]
    vntr_region = config["bam_processing"]["vntr_region_hg38"]
```

### 3.2 VNtyper Strengths

✅ **Simple and explicit** - Clear mapping between assembly and resources
✅ **Config-driven** - No hardcoded paths in code
✅ **Assembly-aware** - Each tool/module can have assembly-specific paths
✅ **User-friendly** - Single `--reference-assembly` flag controls everything

### 3.3 VNtyper Limitations

⚠️ **Not extensible** - Adding new assembly requires config schema changes
⚠️ **No validation** - Doesn't check if paths exist or match assembly
⚠️ **Scattered logic** - If/else blocks throughout codebase
⚠️ **No provenance** - Can't track reference source/version
⚠️ **Hardcoded assemblies** - Only hg19/hg38 supported

---

## 4. Critical Requirements

### 4.1 Functional Requirements

**FR1: Multi-Assembly Support**
- Support hg19, hg38, GRCh37, GRCh38, and future assemblies
- Enable custom assemblies for specialized research
- Support non-human organisms (mouse, etc.)

**FR2: Path Management**
- Centralized mapping: assembly → FASTA paths
- Support multiple references per assembly (with/without ALT, decoy)
- Automatic index discovery (BWA, minimap2, samtools)

**FR3: Validation**
- Verify reference file existence before pipeline execution
- Check index integrity (BWA, minimap2)
- Validate assembly consistency across runs
- Warn about common issues (missing indices, wrong assembly)

**FR4: Chromosome Naming Handling**
- Auto-detect UCSC vs NCBI naming conventions
- Convert coordinates between conventions when needed
- Handle mitochondrial sequence variations

**FR5: Provenance & Reproducibility**
- Record reference source (URL, version)
- Compute reference checksums (MD5/SHA256)
- Track which assembly was used in simulation outputs

### 4.2 Non-Functional Requirements

**NFR1: Backward Compatibility**
- Existing configs must continue working
- Gradual migration path for users
- Legacy format auto-conversion

**NFR2: Simplicity**
- No external dependencies (avoid Refgenie dependency)
- Config-driven (no database required)
- Intuitive for bioinformaticians

**NFR3: Performance**
- Reference validation cached
- Lazy loading (validate only when needed)
- No startup overhead

**NFR4: Extensibility**
- Easy to add new assemblies
- Support organism-specific config
- Plugin architecture for custom validators

---

## 5. Comprehensive Solution Design

### 5.1 Enhanced Config Schema

**New `reference_genomes` Section**:
```json
{
  "reference_genomes": {
    "hg38": {
      "display_name": "Human GRCh38/hg38 (UCSC)",
      "organism": "homo_sapiens",
      "assembly_name": "GRCh38",
      "assembly_alias": ["hg38", "GRCh38.p14"],
      "chromosome_style": "ucsc",  // "ucsc" (chr1) or "ncbi" (1)
      "fasta_path": "${REFERENCE_DIR}/hg38/hg38.fa",
      "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
      "fasta_md5": "TODO",
      "indices": {
        "bwa": "${REFERENCE_DIR}/hg38/bwa/",
        "minimap2": "${REFERENCE_DIR}/hg38/hg38.fa.mmi",
        "samtools": "${REFERENCE_DIR}/hg38/hg38.fa.fai"
      },
      "muc1_vntr_region": "chr1:155188487-155192239",
      "metadata": {
        "release_date": "2013-12-17",
        "source": "UCSC Genome Browser",
        "patch_version": "p14",
        "notes": "Recommended for new projects (2024+)"
      }
    },
    "hg19": {
      "display_name": "Human GRCh37/hg19 (UCSC)",
      "organism": "homo_sapiens",
      "assembly_name": "GRCh37",
      "assembly_alias": ["hg19", "b37"],
      "chromosome_style": "ucsc",
      "fasta_path": "${REFERENCE_DIR}/hg19/hg19.fa",
      "fasta_url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
      "fasta_md5": "TODO",
      "indices": {
        "bwa": "${REFERENCE_DIR}/hg19/bwa/",
        "minimap2": "${REFERENCE_DIR}/hg19/hg19.fa.mmi",
        "samtools": "${REFERENCE_DIR}/hg19/hg19.fa.fai"
      },
      "muc1_vntr_region": "chr1:155160963-155162030",
      "metadata": {
        "release_date": "2009-02-27",
        "source": "UCSC Genome Browser",
        "notes": "Legacy assembly, use hg38 for new projects"
      }
    },
    "GRCh38": {
      "display_name": "Human GRCh38 (NCBI)",
      "organism": "homo_sapiens",
      "assembly_name": "GRCh38",
      "assembly_alias": ["GRCh38.p14"],
      "chromosome_style": "ncbi",  // No "chr" prefix
      "fasta_path": "${REFERENCE_DIR}/GRCh38/GRCh38.fa",
      "fasta_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz",
      "fasta_md5": "TODO",
      "indices": {
        "bwa": "${REFERENCE_DIR}/GRCh38/bwa/",
        "minimap2": "${REFERENCE_DIR}/GRCh38/GRCh38.fa.mmi",
        "samtools": "${REFERENCE_DIR}/GRCh38/GRCh38.fa.fai"
      },
      "muc1_vntr_region": "1:155188487-155192239",  // No "chr" prefix
      "metadata": {
        "release_date": "2013-12-17",
        "source": "NCBI/GRC",
        "patch_version": "p14"
      }
    },
    "mm39": {
      "display_name": "Mouse GRCm39/mm39 (UCSC)",
      "organism": "mus_musculus",
      "assembly_name": "GRCm39",
      "assembly_alias": ["mm39"],
      "chromosome_style": "ucsc",
      "fasta_path": "${REFERENCE_DIR}/mm39/mm39.fa",
      "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
      "indices": {},
      "metadata": {
        "release_date": "2020-06",
        "source": "UCSC Genome Browser",
        "notes": "Recommended for mouse (2024+)"
      }
    }
  },
  "reference_assembly": "hg38",  // Active assembly (backward compatible)
  "constants": {
    // Preserved for backward compatibility
    // Auto-populated from reference_genomes[reference_assembly]
    "hg38": { ... },
    "hg19": { ... }
  }
}
```

**Environment Variable Support**:
- `${REFERENCE_DIR}` → Expands to configured reference directory
- `${HOME}` → User home directory
- `${MUCONEUP_DATA}` → XDG data directory

### 5.2 New Reference Manager Module

**Location**: `muc_one_up/bioinformatics/reference_manager.py`

```python
"""Reference genome management for MucOneUp.

This module provides a centralized interface for managing reference genomes,
including path resolution, validation, and metadata access.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..exceptions import FileOperationError, ValidationError
from ..type_defs import ConfigDict, FilePath
from .reference_validation import validate_reference_genome


@dataclass
class ReferenceGenome:
    """Reference genome metadata and paths.

    Attributes:
        name: Short name (e.g., "hg38", "hg19")
        display_name: Human-readable name
        organism: Organism name (homo_sapiens, mus_musculus)
        assembly_name: Official assembly name (GRCh38, GRCh37)
        assembly_alias: Alternative names
        chromosome_style: "ucsc" (chr1) or "ncbi" (1)
        fasta_path: Path to FASTA file
        fasta_url: Download URL for FASTA
        fasta_md5: MD5 checksum
        indices: Dictionary of index paths
        vntr_region: MUC1 VNTR region coordinates
        metadata: Additional metadata
    """
    name: str
    display_name: str
    organism: str
    assembly_name: str
    assembly_alias: list[str]
    chromosome_style: str
    fasta_path: Path
    fasta_url: str | None
    fasta_md5: str | None
    indices: dict[str, str]
    vntr_region: str | None
    metadata: dict[str, Any]

    @property
    def is_ucsc_style(self) -> bool:
        """Check if chromosome naming uses UCSC convention (chr1)."""
        return self.chromosome_style.lower() == "ucsc"

    @property
    def is_ncbi_style(self) -> bool:
        """Check if chromosome naming uses NCBI convention (1)."""
        return self.chromosome_style.lower() == "ncbi"


class ReferenceManager:
    """Centralized reference genome management.

    Provides access to configured reference genomes, validates paths,
    and manages reference-related resources.
    """

    def __init__(self, config: ConfigDict):
        """Initialize reference manager.

        Args:
            config: MucOneUp configuration dictionary
        """
        self.config = config
        self._references: dict[str, ReferenceGenome] = {}
        self._load_references()

    def _expand_path(self, path_template: str) -> Path:
        """Expand environment variables in path template.

        Args:
            path_template: Path with potential variables (${VAR})

        Returns:
            Expanded absolute path
        """
        import os
        import re

        # Replace ${VAR} patterns
        def replace_var(match: re.Match) -> str:
            var_name = match.group(1)
            if var_name == "REFERENCE_DIR":
                return self.config.get("reference_directory", "./reference")
            return os.environ.get(var_name, match.group(0))

        expanded = re.sub(r'\$\{(\w+)\}', replace_var, path_template)
        return Path(expanded).expanduser().resolve()

    def _load_references(self) -> None:
        """Load reference genome configurations."""
        ref_genomes = self.config.get("reference_genomes", {})

        for ref_name, ref_config in ref_genomes.items():
            fasta_path = self._expand_path(ref_config["fasta_path"])

            # Expand index paths
            indices = {}
            for idx_type, idx_path in ref_config.get("indices", {}).items():
                indices[idx_type] = str(self._expand_path(idx_path))

            self._references[ref_name] = ReferenceGenome(
                name=ref_name,
                display_name=ref_config.get("display_name", ref_name),
                organism=ref_config.get("organism", "unknown"),
                assembly_name=ref_config.get("assembly_name", ref_name),
                assembly_alias=ref_config.get("assembly_alias", []),
                chromosome_style=ref_config.get("chromosome_style", "ucsc"),
                fasta_path=fasta_path,
                fasta_url=ref_config.get("fasta_url"),
                fasta_md5=ref_config.get("fasta_md5"),
                indices=indices,
                vntr_region=ref_config.get("muc1_vntr_region"),
                metadata=ref_config.get("metadata", {})
            )

    def get_reference(self, assembly: str | None = None) -> ReferenceGenome:
        """Get reference genome by name.

        Args:
            assembly: Assembly name. If None, uses config["reference_assembly"]

        Returns:
            ReferenceGenome object

        Raises:
            ValidationError: If assembly not found
        """
        if assembly is None:
            assembly = self.config.get("reference_assembly", "hg38")

        # Check direct name
        if assembly in self._references:
            return self._references[assembly]

        # Check aliases
        for ref in self._references.values():
            if assembly in ref.assembly_alias:
                return ref

        available = list(self._references.keys())
        raise ValidationError(
            f"Reference assembly '{assembly}' not found. "
            f"Available: {', '.join(available)}"
        )

    def get_fasta_path(self, assembly: str | None = None) -> Path:
        """Get FASTA file path for assembly.

        Args:
            assembly: Assembly name

        Returns:
            Path to FASTA file

        Raises:
            ValidationError: If assembly not found
            FileOperationError: If FASTA file missing
        """
        ref = self.get_reference(assembly)

        if not ref.fasta_path.exists():
            raise FileOperationError(
                f"Reference FASTA not found: {ref.fasta_path}\n"
                f"Download from: {ref.fasta_url}"
            )

        return ref.fasta_path

    def get_vntr_region(self, assembly: str | None = None) -> str:
        """Get MUC1 VNTR region coordinates.

        Args:
            assembly: Assembly name

        Returns:
            Genomic region string (e.g., "chr1:155188487-155192239")

        Raises:
            ValidationError: If region not defined
        """
        ref = self.get_reference(assembly)

        if ref.vntr_region is None:
            raise ValidationError(
                f"MUC1 VNTR region not defined for assembly: {ref.name}"
            )

        return ref.vntr_region

    def validate_reference(
        self,
        assembly: str | None = None,
        aligner: str = "bwa",
        check_indices: bool = True
    ) -> list[str]:
        """Validate reference genome and indices.

        Args:
            assembly: Assembly name
            aligner: Aligner to check ("bwa" or "minimap2")
            check_indices: Whether to validate indices

        Returns:
            List of warnings (empty if all OK)

        Raises:
            ValidationError: If assembly not found
            FileOperationError: If FASTA missing
        """
        ref = self.get_reference(assembly)
        fasta_path = self.get_fasta_path(assembly)

        warnings: list[str] = []

        if check_indices:
            warnings = validate_reference_genome(fasta_path, aligner=aligner)

        return warnings

    def list_assemblies(
        self,
        organism: str | None = None
    ) -> list[str]:
        """List available assemblies.

        Args:
            organism: Filter by organism (e.g., "homo_sapiens")

        Returns:
            List of assembly names
        """
        assemblies = []
        for ref in self._references.values():
            if organism is None or ref.organism == organism:
                assemblies.append(ref.name)
        return assemblies

    def convert_chromosome_name(
        self,
        chrom: str,
        from_style: str,
        to_style: str
    ) -> str:
        """Convert chromosome name between UCSC and NCBI conventions.

        Args:
            chrom: Chromosome name
            from_style: Source style ("ucsc" or "ncbi")
            to_style: Target style ("ucsc" or "ncbi")

        Returns:
            Converted chromosome name

        Examples:
            >>> convert_chromosome_name("chr1", "ucsc", "ncbi")
            "1"
            >>> convert_chromosome_name("MT", "ncbi", "ucsc")
            "chrM"
        """
        if from_style == to_style:
            return chrom

        if from_style == "ucsc" and to_style == "ncbi":
            # Remove "chr" prefix
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            # Handle mitochondrial
            if chrom == "M":
                return "MT"
            return chrom

        elif from_style == "ncbi" and to_style == "ucsc":
            # Add "chr" prefix
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            # Handle mitochondrial
            if chrom == "chrMT":
                return "chrM"
            return chrom

        return chrom


# Convenience function for backward compatibility
def get_reference_path_for_assembly(
    config: ConfigDict,
    assembly: str | None = None
) -> Path:
    """Get reference FASTA path for assembly.

    This is a convenience function that maintains backward compatibility
    with the original issue #28 proposal.

    Args:
        config: Configuration dictionary
        assembly: Assembly name

    Returns:
        Path to reference FASTA

    Raises:
        ValidationError: If assembly not configured
        FileOperationError: If reference missing
    """
    manager = ReferenceManager(config)
    return manager.get_fasta_path(assembly)
```

### 5.3 Integration Points

**Pipeline Integration**:
```python
# muc_one_up/read_simulator/ont_pipeline.py

from ..bioinformatics.reference_manager import ReferenceManager

def simulate_ont_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None
) -> str:
    # If human_reference not provided, get from reference manager
    if human_reference is None:
        try:
            ref_manager = ReferenceManager(config)
            assembly = config.get("reference_assembly", "hg38")
            human_reference = str(ref_manager.get_fasta_path(assembly))

            # Validate indices
            warnings = ref_manager.validate_reference(assembly, aligner="minimap2")
            for warning in warnings:
                logging.warning(warning)

            logging.info(f"Using reference: {human_reference} ({assembly})")
        except (ValidationError, FileOperationError) as e:
            logging.warning(f"Reference not configured: {e}")
            logging.warning("Will align to simulated reference instead")
            human_reference = input_fa

    # Rest of pipeline...
```

**CLI Integration**:
```python
# muc_one_up/cli/click_main.py

@click.command()
@click.option("--reference-assembly",
              type=str,
              help="Reference assembly to use (e.g., hg38, hg19, GRCh38)")
def simulate(reference_assembly, ...):
    # Load config
    config = load_config(config_path)

    # Override reference assembly if provided
    if reference_assembly:
        config["reference_assembly"] = reference_assembly

    # Validate reference before starting
    ref_manager = ReferenceManager(config)
    ref = ref_manager.get_reference()

    click.echo(f"Using reference: {ref.display_name}")

    warnings = ref_manager.validate_reference()
    if warnings:
        click.echo("⚠️  Reference warnings:")
        for warning in warnings:
            click.echo(f"  - {warning}")

    # Continue with simulation...
```

### 5.4 Enhanced Helper Script

**`helpers/manage_references.py`** (New):
```python
#!/usr/bin/env python3
"""
Reference genome management utility for MucOneUp.

Provides commands for downloading, validating, and managing reference genomes
and their associated indices.
"""

import click
from pathlib import Path

from muc_one_up.config import load_config
from muc_one_up.bioinformatics.reference_manager import ReferenceManager


@click.group()
def cli():
    """MucOneUp reference genome management."""
    pass


@cli.command()
@click.option("--config", type=Path, default="config.json",
              help="Path to config file")
def list(config):
    """List available reference assemblies."""
    cfg = load_config(str(config))
    manager = ReferenceManager(cfg)

    click.echo("Available reference assemblies:\n")

    for name in manager.list_assemblies():
        ref = manager.get_reference(name)
        exists = "✓" if ref.fasta_path.exists() else "✗"
        click.echo(f"  [{exists}] {ref.name:10s} - {ref.display_name}")
        click.echo(f"      {ref.fasta_path}")
        if ref.vntr_region:
            click.echo(f"      VNTR: {ref.vntr_region}")
        click.echo()


@cli.command()
@click.option("--config", type=Path, default="config.json")
@click.argument("assembly")
@click.option("--aligner", type=click.Choice(["bwa", "minimap2"]), default="bwa")
def validate(config, assembly, aligner):
    """Validate reference genome and indices."""
    cfg = load_config(str(config))
    manager = ReferenceManager(cfg)

    ref = manager.get_reference(assembly)
    click.echo(f"Validating {ref.display_name}...")

    try:
        warnings = manager.validate_reference(assembly, aligner=aligner)

        if warnings:
            click.echo("\n⚠️  Warnings:")
            for warning in warnings:
                click.echo(f"  - {warning}")
        else:
            click.echo("✓ All checks passed")
    except Exception as e:
        click.echo(f"✗ Validation failed: {e}", err=True)
        raise SystemExit(1)


@cli.command()
@click.option("--config", type=Path, default="config.json")
@click.argument("assembly")
def download(config, assembly):
    """Download reference genome."""
    cfg = load_config(str(config))
    manager = ReferenceManager(cfg)

    ref = manager.get_reference(assembly)

    if ref.fasta_url is None:
        click.echo(f"No download URL configured for {assembly}", err=True)
        raise SystemExit(1)

    click.echo(f"Downloading {ref.display_name}...")
    click.echo(f"Source: {ref.fasta_url}")

    # Implementation using download_references.py logic
    # ...


if __name__ == "__main__":
    cli()
```

---

## 6. Implementation Plan

### Phase 1: Foundation (Week 1)

**Tasks**:
1. ✅ Create `reference_manager.py` module
2. ✅ Implement `ReferenceGenome` dataclass
3. ✅ Implement `ReferenceManager` class
4. ✅ Add path expansion (${REFERENCE_DIR})
5. ✅ Add unit tests for reference manager

**Deliverables**:
- `muc_one_up/bioinformatics/reference_manager.py`
- `tests/bioinformatics/test_reference_manager.py`

### Phase 2: Config Schema (Week 1)

**Tasks**:
1. ✅ Update `CONFIG_SCHEMA` in `config.py`
2. ✅ Add `reference_genomes` section validation
3. ✅ Implement backward compatibility conversion
4. ✅ Create example `config.json` with all assemblies
5. ✅ Update schema documentation

**Deliverables**:
- Updated `muc_one_up/config.py`
- Updated `config.json` (example)
- Migration guide

### Phase 3: Pipeline Integration (Week 2)

**Tasks**:
1. ✅ Update `ont_pipeline.py` to use ReferenceManager
2. ✅ Update `pipeline.py` (Illumina) to use ReferenceManager
3. ✅ Add reference validation to CLI
4. ✅ Update `--reference-assembly` CLI option
5. ✅ Add integration tests

**Deliverables**:
- Updated pipelines
- Updated CLI commands
- Integration tests

### Phase 4: Helper Tools (Week 2)

**Tasks**:
1. ✅ Create `manage_references.py` CLI tool
2. ✅ Update `download_references.py` to use new schema
3. ✅ Add reference listing command
4. ✅ Add validation command
5. ✅ Add download command

**Deliverables**:
- `helpers/manage_references.py`
- Updated `helpers/download_references.py`

### Phase 5: Documentation (Week 3)

**Tasks**:
1. ✅ Update CLAUDE.md with reference management
2. ✅ Update README.md with examples
3. ✅ Create reference configuration guide
4. ✅ Document chromosome naming conventions
5. ✅ Create migration guide for existing users

**Deliverables**:
- Updated documentation
- Migration guide
- Configuration examples

### Phase 6: Testing & Validation (Week 3)

**Tasks**:
1. ✅ Full integration testing with all assemblies
2. ✅ Test backward compatibility
3. ✅ Test chromosome name conversion
4. ✅ Performance benchmarking
5. ✅ User acceptance testing

**Deliverables**:
- Test coverage report
- Performance benchmarks
- UAT results

---

## 7. Migration Strategy

### 7.1 Backward Compatibility

**Existing Config Format** (still works):
```json
{
  "constants": {
    "hg38": { "left": "...", "right": "...", "vntr_region": "..." },
    "hg19": { "left": "...", "right": "...", "vntr_region": "..." }
  },
  "reference_assembly": "hg38"
}
```

**Auto-Migration**:
- `load_config()` detects old format
- Automatically creates minimal `reference_genomes` section
- Logs migration notice
- Users encouraged to upgrade to new format

### 7.2 Migration Steps for Users

1. **Check current assembly**:
   ```bash
   muconeup --config config.json reference list
   ```

2. **Validate current reference**:
   ```bash
   muconeup --config config.json reference validate hg38
   ```

3. **Update config** (optional):
   - Add `reference_genomes` section from template
   - Set `REFERENCE_DIR` environment variable
   - Validate with `muconeup --config config.json validate`

4. **Download missing references**:
   ```bash
   muconeup --config config.json reference download hg38
   ```

---

## 8. Testing Strategy

### 8.1 Unit Tests

**Test Coverage**:
- ✅ `ReferenceManager` initialization
- ✅ Reference loading from config
- ✅ Path expansion (environment variables)
- ✅ Reference retrieval (by name, by alias)
- ✅ Error handling (missing assembly, missing file)
- ✅ Chromosome name conversion
- ✅ VNTR region retrieval

### 8.2 Integration Tests

**Test Scenarios**:
- ✅ Simulate with hg38 reference
- ✅ Simulate with hg19 reference
- ✅ Switch assembly mid-workflow (should fail with warning)
- ✅ Run with custom assembly
- ✅ Run with missing indices (should warn)

### 8.3 Validation Tests

**Test Matrix**:

| Assembly | Chr Style | Exists | Indices | Expected Result |
|----------|-----------|--------|---------|-----------------|
| hg38     | UCSC      | ✓      | ✓       | Success         |
| hg38     | UCSC      | ✓      | ✗       | Warning         |
| hg38     | UCSC      | ✗      | -       | Error           |
| hg19     | UCSC      | ✓      | ✓       | Success         |
| GRCh38   | NCBI      | ✓      | ✓       | Success         |
| custom   | UCSC      | ✓      | ✗       | Warning         |

---

## 9. Documentation Requirements

### 9.1 User-Facing Documentation

**README.md**:
- Quick start with default assembly (hg38)
- How to switch assemblies
- How to add custom assemblies
- Reference download instructions

**Configuration Guide**:
- `reference_genomes` schema
- Environment variable usage
- Chromosome naming conventions
- Index management

**Migration Guide**:
- Upgrading from old config format
- Migrating existing analyses
- Troubleshooting common issues

### 9.2 Developer Documentation

**CLAUDE.md**:
- Reference management architecture
- How to extend for new organisms
- Testing reference handling
- API reference for `ReferenceManager`

**Code Documentation**:
- Docstrings for all public functions
- Type hints throughout
- Examples in docstrings
- Architecture diagrams

---

## Appendix A: Reference Genome Metadata

### A.1 Human Assemblies

| Assembly | Release | Chromosome Style | MUC1 VNTR Region | Notes |
|----------|---------|------------------|------------------|-------|
| hg38 (GRCh38) | 2013-12 | UCSC (chr1) | chr1:155188487-155192239 | **Recommended** |
| hg19 (GRCh37) | 2009-02 | UCSC (chr1) | chr1:155160963-155162030 | Legacy |
| GRCh38.p14 | 2022-02 | NCBI (1) | 1:155188487-155192239 | Latest patch |
| GRCh37.p13 | 2013-06 | NCBI (1) | 1:155160963-155162030 | Final GRCh37 |

### A.2 Mouse Assemblies

| Assembly | Release | Chromosome Style | Notes |
|----------|---------|------------------|-------|
| mm39 (GRCm39) | 2020-06 | UCSC (chr1) | **Recommended** |
| mm10 (GRCm38) | 2011-12 | UCSC (chr1) | Legacy |

---

## Appendix B: Comparison Matrix

| Feature | Current | VNtyper | Refgenie | Proposed |
|---------|---------|---------|----------|----------|
| **Multi-assembly** | ✗ | ✗ | ✓ | ✓ |
| **Path management** | ✗ | ✓ | ✓ | ✓ |
| **Validation** | Partial | ✗ | ✓ | ✓ |
| **Chr naming** | ✗ | ✗ | ✓ | ✓ |
| **Provenance** | ✗ | ✗ | ✓ | ✓ |
| **Backward compat** | - | - | - | ✓ |
| **Simplicity** | ✓ | ✓ | ✗ | ✓ |
| **No dependencies** | ✓ | ✓ | ✗ | ✓ |
| **Extensibility** | ✗ | ✗ | ✓ | ✓ |

---

## Appendix C: Example Workflows

### C.1 Standard Workflow (hg38)

```bash
# 1. Validate reference
muconeup --config config.json reference validate hg38

# 2. Run simulation
muconeup --config config.json simulate \
  --reference-assembly hg38 \
  --out-base sample_001

# 3. Generate reads
muconeup --config config.json reads ont \
  sample_001.fa \
  --seed 42
```

### C.2 Custom Assembly Workflow

```json
// Add to config.json
{
  "reference_genomes": {
    "custom_hg38": {
      "display_name": "Custom hg38 with decoy",
      "organism": "homo_sapiens",
      "assembly_name": "GRCh38",
      "chromosome_style": "ucsc",
      "fasta_path": "/data/custom/hg38_with_decoy.fa",
      "muc1_vntr_region": "chr1:155188487-155192239"
    }
  }
}
```

```bash
# Use custom assembly
muconeup --config config.json simulate \
  --reference-assembly custom_hg38 \
  --out-base sample_custom
```

---

## Conclusion

This comprehensive solution provides:

✅ **Robust multi-assembly support** - Easily add new assemblies
✅ **Industry best practices** - Follows GATK/Refgenie philosophies
✅ **Backward compatibility** - Existing configs continue working
✅ **Simplicity** - No external dependencies, config-driven
✅ **Extensibility** - Ready for future organisms and assemblies
✅ **Validation** - Catches errors before pipeline execution
✅ **Provenance** - Track reference sources and versions

The design balances enterprise-grade reference management with MucOneUp's philosophy of simplicity and ease of use. By learning from VNtyper's pragmatic approach and Refgenie's comprehensive asset management, we achieve the best of both worlds.

**Next Steps**:
1. Review and approve design
2. Begin Phase 1 implementation
3. Iterate based on feedback
4. Deploy with thorough testing

---

**End of Assessment**

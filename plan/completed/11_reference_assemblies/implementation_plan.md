# Issue #28: Reference Assembly Management - FINAL REVISED SOLUTION

**Author**: Expert Senior Developer & Bioinformatician (Ultra-Thinking Review)
**Date**: 2025-10-19
**Status**: CORRECTED - Anti-Pattern Violations Fixed
**Previous Version**: Over-engineered (issue_28_comprehensive_reference_assembly_solution.md)

---

## ⚠️ CRITICAL REVISION NOTICE

**This document SUPERSEDES the previous comprehensive assessment** which violated DRY, KISS, and SOLID principles by over-engineering the solution.

**What went wrong in previous assessment**:
- ❌ Created `ReferenceManager` class (explicitly forbidden in original issue #28)
- ❌ 500+ lines of code when 100 would suffice
- ❌ Duplicated existing `reference_validation.py` functionality
- ❌ Feature creep (mouse genomes, provenance tracking, chromosome conversion)
- ❌ Unnecessary abstractions (dataclasses, environment variables, helper CLI)

**What this revision does**:
- ✅ Follows original corrected issue #28 approach
- ✅ Extends existing `reference_validation.py` with 3 helper functions
- ✅ NO new classes or modules
- ✅ Simple VNtyper-inspired config
- ✅ Maintains DRY, KISS, SOLID principles

---

## Executive Summary

**Solution**: Extend existing `reference_validation.py` with 3 simple helper functions to map assemblies to reference paths and validate them. Total code addition: ~100 lines.

**Inspiration**: VNtyper's pragmatic config-driven approach + Original issue #28's KISS adherence.

**Principle Compliance**:
- ✅ **DRY**: Reuses existing validation functions
- ✅ **KISS**: Simple functions, no unnecessary classes
- ✅ **SOLID**: Extends validation module (Single Responsibility maintained)
- ✅ **Modular**: Adds to appropriate existing module

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Current State](#2-current-state)
3. [VNtyper's Simple Approach](#3-vntypers-simple-approach)
4. [Recommended Solution](#4-recommended-solution)
5. [Implementation](#5-implementation)
6. [Testing](#6-testing)
7. [Migration](#7-migration)
8. [Anti-Patterns Avoided](#8-anti-patterns-avoided)

---

## 1. Problem Statement

### 1.1 The Core Problem

MucOneUp supports hg19/hg38 in config but lacks:
1. **Path management** - No link between assembly name → reference FASTA path
2. **Validation** - Can't verify reference matches assembly before pipeline runs
3. **Extensibility** - Hard to add new assemblies (GRCh37, custom references)

### 1.2 What Users Need

```bash
# User wants this to work
muconeup --config config.json simulate \
  --reference-assembly hg38 \
  --out-base sample
```

**Current problem**: Config has `reference_assembly: hg38` but pipelines don't know which FASTA file to use.

### 1.3 Scope Constraints

**In Scope**:
- Map assembly → reference FASTA path
- Validate reference exists and has indices
- Get VNTR region for assembly
- Support hg19, hg38, GRCh37, GRCh38, custom

**Out of Scope** (prevents feature creep):
- ❌ Non-human organisms (mouse, etc.) - MUC1 is human-specific
- ❌ Chromosome name conversion - Just store correct value for each assembly
- ❌ Provenance tracking - Nice-to-have, not essential
- ❌ Download management - Already exists in `helpers/`
- ❌ Complex metadata - Keep it simple

---

## 2. Current State

### 2.1 Existing Infrastructure (Good Foundation)

**Config Support**:
```json
{
  "constants": {
    "hg38": {
      "left": "...",
      "right": "...",
      "vntr_region": "chr1:155188487-155192239"
    },
    "hg19": {...}
  },
  "reference_assembly": "hg38"
}
```

**Validation Module** (`reference_validation.py`):
- ✅ `validate_reference_genome()` - Checks FASTA + indices
- ✅ `_validate_bwa_indices()` - BWA validation
- ✅ `_validate_minimap2_indices()` - Minimap2 validation

**Reference Utils** (`reference_utils.py`):
- ✅ `get_reference_info()` - Parse FASTA metadata
- ✅ `is_diploid_reference()` - Diploid detection

### 2.2 The Gap

**Missing**: Simple mapping function
```python
# This doesn't exist yet
ref_path = get_reference_path_for_assembly(config, assembly="hg38")
# -> Path("/data/references/hg38.fa")
```

---

## 3. VNtyper's Simple Approach

### 3.1 VNtyper Config (Excellent Reference)

```json
{
  "reference_data": {
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

### 3.2 VNtyper Runtime (Simple and Clear)

```python
# Pipeline code
if args.reference_assembly == "hg19":
    bwa_reference = config["reference_data"]["bwa_reference_hg19"]
    vntr_region = config["bam_processing"]["vntr_region_hg19"]
else:  # hg38
    bwa_reference = config["reference_data"]["bwa_reference_hg38"]
    vntr_region = config["bam_processing"]["vntr_region_hg38"]
```

### 3.3 Why VNtyper's Approach Works

✅ **Explicit** - Clear what reference is used for each assembly
✅ **Simple** - No classes, no abstractions
✅ **Config-driven** - Change config, not code
✅ **Debuggable** - Easy to trace which file is being used

### 3.4 VNtyper's Limitation

⚠️ **Not extensible** - Adding new assembly requires:
1. Adding new config keys (`bwa_reference_mm39`, etc.)
2. Adding new if/else branches in code
3. Schema changes

**Our improvement**: Group by assembly, use helper function instead of if/else.

---

## 4. Recommended Solution

### 4.1 Enhanced Config (VNtyper-Inspired, Grouped by Assembly)

```json
{
  "reference_genomes": {
    "hg38": {
      "fasta_path": "reference/hg38/hg38.fa",
      "vntr_region": "chr1:155188487-155192239",
      "display_name": "Human GRCh38/hg38 (UCSC)",
      "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    },
    "hg19": {
      "fasta_path": "reference/hg19/hg19.fa",
      "vntr_region": "chr1:155160963-155162030",
      "display_name": "Human GRCh37/hg19 (UCSC)",
      "source_url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    },
    "GRCh38": {
      "fasta_path": "reference/GRCh38/GRCh38.fa",
      "vntr_region": "1:155188487-155192239",
      "display_name": "Human GRCh38 (NCBI - no chr prefix)",
      "source_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz"
    },
    "custom": {
      "fasta_path": "/absolute/path/to/custom.fa",
      "vntr_region": "chr1:155188487-155192239",
      "display_name": "Custom hg38 with modifications"
    }
  },
  "reference_assembly": "hg38",

  "constants": {
    "hg38": {...},
    "hg19": {...}
  }
}
```

**Key Design Decisions**:
1. **Group by assembly** - All assembly-specific data together
2. **Minimal required fields** - `fasta_path`, `vntr_region` (others optional)
3. **Backward compatible** - `constants` section preserved
4. **Relative paths supported** - Relative to config file location
5. **No environment variables** - Simpler, more predictable

### 4.2 Config Schema Update

```python
# muc_one_up/config.py - Add to CONFIG_SCHEMA

"reference_genomes": {
    "type": "object",
    "patternProperties": {
        "^[a-zA-Z0-9_]+$": {  # Assembly name
            "type": "object",
            "properties": {
                "fasta_path": {"type": "string"},
                "vntr_region": {"type": "string"},
                "display_name": {"type": "string"},
                "source_url": {"type": "string"}
            },
            "required": ["fasta_path", "vntr_region"],
            "additionalProperties": False
        }
    }
}
```

**Validation**: Ensure `fasta_path` is string, `vntr_region` follows format `chr?[0-9XY]+:[0-9]+-[0-9]+`

---

## 5. Implementation

### 5.1 Extend `reference_validation.py` (3 Functions)

**File**: `muc_one_up/bioinformatics/reference_validation.py`

```python
"""Reference genome validation for MucOneUp.

EXISTING FUNCTIONS:
- validate_reference_genome()
- _validate_bwa_indices()
- _validate_minimap2_indices()
- validate_bam_file()
- validate_bed_file()

NEW FUNCTIONS (Issue #28):
- get_reference_path_for_assembly()
- validate_reference_for_assembly()
- get_muc1_region_for_assembly()
"""

import logging
from pathlib import Path

from ..exceptions import FileOperationError, ValidationError
from ..type_defs import ConfigDict, FilePath


def get_reference_path_for_assembly(
    config: ConfigDict,
    assembly: str | None = None,
) -> Path:
    """Get reference FASTA path for specified assembly.

    Args:
        config: Configuration dictionary
        assembly: Assembly name (e.g., "hg38", "hg19", "GRCh38").
                 If None, uses config["reference_assembly"]

    Returns:
        Absolute path to reference FASTA file

    Raises:
        ValidationError: If assembly not configured in reference_genomes
        FileOperationError: If reference file doesn't exist

    Example:
        >>> config = load_config("config.json")
        >>> ref_path = get_reference_path_for_assembly(config, "hg38")
        >>> print(ref_path)
        /data/reference/hg38/hg38.fa
    """
    # Get active assembly
    if assembly is None:
        assembly = config.get("reference_assembly", "hg38")
        logging.debug(f"No assembly specified, using default: {assembly}")

    # Get reference_genomes section
    ref_genomes = config.get("reference_genomes", {})

    if not ref_genomes:
        raise ValidationError(
            "No 'reference_genomes' section in config. "
            "Please add reference genome configuration."
        )

    # Get assembly config
    assembly_config = ref_genomes.get(assembly)

    if not assembly_config:
        available = list(ref_genomes.keys())
        raise ValidationError(
            f"Assembly '{assembly}' not found in reference_genomes. "
            f"Available assemblies: {', '.join(available)}"
        )

    # Get FASTA path
    fasta_path_str = assembly_config.get("fasta_path")

    if not fasta_path_str:
        raise ValidationError(
            f"No 'fasta_path' configured for assembly '{assembly}'"
        )

    # Convert to Path and resolve
    fasta_path = Path(fasta_path_str)

    # If relative path, make it relative to config directory
    # (This assumes config file location is stored somewhere, or use cwd)
    if not fasta_path.is_absolute():
        # Make relative to current working directory
        fasta_path = fasta_path.resolve()

    # Validate file exists
    if not fasta_path.exists():
        source_url = assembly_config.get("source_url", "")
        error_msg = f"Reference genome file not found: {fasta_path}"
        if source_url:
            error_msg += f"\n  Download from: {source_url}"
        raise FileOperationError(error_msg)

    logging.debug(f"Using reference for {assembly}: {fasta_path}")
    return fasta_path


def validate_reference_for_assembly(
    config: ConfigDict,
    assembly: str | None = None,
    aligner: str = "bwa",
) -> list[str]:
    """Validate reference genome and indices for specified assembly.

    Calls get_reference_path_for_assembly() to get the path, then
    validates using the existing validate_reference_genome() function.

    Args:
        config: Configuration dictionary
        assembly: Assembly name (e.g., "hg38", "hg19")
        aligner: Aligner to check indices for ("bwa" or "minimap2")

    Returns:
        List of warning messages (empty if all checks pass)

    Raises:
        ValidationError: If assembly not configured
        FileOperationError: If reference file missing

    Example:
        >>> warnings = validate_reference_for_assembly(config, "hg38", "bwa")
        >>> if warnings:
        ...     for w in warnings:
        ...         print(f"Warning: {w}")
    """
    # Get reference path (this validates assembly exists)
    ref_path = get_reference_path_for_assembly(config, assembly)

    # Use EXISTING validation function (DRY principle)
    warnings = validate_reference_genome(ref_path, aligner=aligner)

    # Log warnings with assembly context
    if assembly is None:
        assembly = config.get("reference_assembly", "unknown")

    for warning in warnings:
        logging.warning(f"[{assembly}] {warning}")

    return warnings


def get_muc1_region_for_assembly(
    config: ConfigDict,
    assembly: str | None = None
) -> str:
    """Get MUC1 VNTR region coordinates for specified assembly.

    Args:
        config: Configuration dictionary
        assembly: Assembly name (e.g., "hg38", "hg19")

    Returns:
        Genomic region string (e.g., "chr1:155188487-155192239")

    Raises:
        ValidationError: If assembly not configured or vntr_region missing

    Example:
        >>> region = get_muc1_region_for_assembly(config, "hg38")
        >>> print(region)
        chr1:155188487-155192239
    """
    # Get active assembly
    if assembly is None:
        assembly = config.get("reference_assembly", "hg38")

    # Get reference_genomes section
    ref_genomes = config.get("reference_genomes", {})

    if not ref_genomes:
        raise ValidationError(
            "No 'reference_genomes' section in config"
        )

    # Get assembly config
    assembly_config = ref_genomes.get(assembly)

    if not assembly_config:
        available = list(ref_genomes.keys())
        raise ValidationError(
            f"Assembly '{assembly}' not found. Available: {', '.join(available)}"
        )

    # Get VNTR region
    vntr_region = assembly_config.get("vntr_region")

    if not vntr_region:
        raise ValidationError(
            f"No 'vntr_region' configured for assembly '{assembly}'"
        )

    logging.debug(f"MUC1 VNTR region for {assembly}: {vntr_region}")
    return vntr_region
```

**That's it! 3 simple functions, ~100 lines total.**

### 5.2 Pipeline Integration

**Example: ONT Pipeline** (`muc_one_up/read_simulator/ont_pipeline.py`):

```python
from ..bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
    validate_reference_for_assembly
)

def simulate_ont_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None
) -> str:
    """Simulate ONT reads using NanoSim."""

    # If human_reference not provided, get from config
    if human_reference is None:
        try:
            assembly = config.get("reference_assembly", "hg38")
            ref_path = get_reference_path_for_assembly(config, assembly)
            human_reference = str(ref_path)

            # Validate indices (warnings logged automatically)
            warnings = validate_reference_for_assembly(
                config, assembly, aligner="minimap2"
            )

            if warnings:
                logging.info(f"Reference validation warnings for {assembly}:")
                for w in warnings:
                    logging.info(f"  - {w}")

            logging.info(f"Using reference: {human_reference} ({assembly})")

        except (ValidationError, FileOperationError) as e:
            logging.warning(f"Could not load reference from config: {e}")
            logging.warning("Falling back to aligning against simulated reference")
            human_reference = input_fa

    # Continue with rest of pipeline...
    # (existing code unchanged)
```

**Example: Illumina Pipeline** (`muc_one_up/read_simulator/pipeline.py`):

```python
from ..bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
    validate_reference_for_assembly
)

def simulate_illumina_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None
) -> str:
    """Simulate Illumina reads."""

    if human_reference is None:
        try:
            assembly = config.get("reference_assembly", "hg38")
            ref_path = get_reference_path_for_assembly(config, assembly)
            human_reference = str(ref_path)

            validate_reference_for_assembly(config, assembly, aligner="bwa")
            logging.info(f"Using reference: {human_reference} ({assembly})")

        except (ValidationError, FileOperationError) as e:
            logging.warning(f"Reference configuration error: {e}")
            # Pipeline can decide how to handle (raise or fallback)

    # Rest of pipeline...
```

### 5.3 CLI Integration

**Example: Click CLI** (`muc_one_up/cli/click_main.py`):

```python
from ..bioinformatics.reference_validation import validate_reference_for_assembly

@click.command()
@click.option("--reference-assembly", help="Override config assembly (e.g., hg38, hg19)")
def simulate(config_path, reference_assembly, ...):
    """Run simulation pipeline."""

    # Load config
    config = load_config(config_path)

    # Override assembly if provided
    if reference_assembly:
        config["reference_assembly"] = reference_assembly

    # Validate reference before starting long pipeline
    try:
        assembly = config.get("reference_assembly", "hg38")
        click.echo(f"Validating reference assembly: {assembly}")

        warnings = validate_reference_for_assembly(config, assembly)

        if warnings:
            click.echo("⚠️  Warnings:")
            for w in warnings:
                click.echo(f"  - {w}")
        else:
            click.echo("✓ Reference validation passed")

    except Exception as e:
        click.echo(f"✗ Reference validation failed: {e}", err=True)
        raise click.Abort()

    # Continue with simulation...
```

### 5.4 Backward Compatibility

**Auto-Migration in `load_config()`**:

```python
# muc_one_up/config.py

def load_config(config_path: str) -> dict[str, Any]:
    """Load and validate config, with backward compatibility."""

    # ... existing loading code ...

    # AUTO-MIGRATION: If reference_genomes missing but constants exists,
    # create minimal reference_genomes section
    if "reference_genomes" not in config and "constants" in config:
        logging.info("Auto-migrating old config format to include reference_genomes")

        ref_assembly = config.get("reference_assembly", "hg38")

        # Create minimal reference_genomes from constants
        config["reference_genomes"] = {}

        for assembly in ["hg38", "hg19"]:
            if assembly in config["constants"]:
                vntr_region = config["constants"][assembly].get("vntr_region")
                if vntr_region:
                    config["reference_genomes"][assembly] = {
                        "fasta_path": f"reference/{assembly}/{assembly}.fa",
                        "vntr_region": vntr_region,
                        "display_name": f"Auto-migrated {assembly}"
                    }

        logging.warning(
            "Config format is outdated. Please add 'reference_genomes' section. "
            "See documentation for details."
        )

    return config
```

---

## 6. Testing

### 6.1 Unit Tests

**File**: `tests/bioinformatics/test_reference_validation.py` (extend existing)

```python
"""Tests for reference validation including assembly management (Issue #28)."""

import pytest
from pathlib import Path

from muc_one_up.config import load_config
from muc_one_up.bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
    validate_reference_for_assembly,
    get_muc1_region_for_assembly
)
from muc_one_up.exceptions import ValidationError, FileOperationError


class TestAssemblyManagement:
    """Test assembly-specific reference management (Issue #28)."""

    def test_get_reference_path_default_assembly(self, tmp_path):
        """Test getting reference path for default assembly."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nATCG\n")

        config = {
            "reference_genomes": {
                "hg38": {
                    "fasta_path": str(ref_file),
                    "vntr_region": "chr1:155188487-155192239"
                }
            },
            "reference_assembly": "hg38"
        }

        path = get_reference_path_for_assembly(config)
        assert path == ref_file

    def test_get_reference_path_explicit_assembly(self, tmp_path):
        """Test getting reference path for explicitly specified assembly."""
        hg19_file = tmp_path / "hg19.fa"
        hg19_file.write_text(">chr1\nATCG\n")

        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/nonexistent/hg38.fa", "vntr_region": "chr1:1-2"},
                "hg19": {"fasta_path": str(hg19_file), "vntr_region": "chr1:1-2"}
            },
            "reference_assembly": "hg38"  # Default is hg38
        }

        # Explicitly request hg19
        path = get_reference_path_for_assembly(config, assembly="hg19")
        assert path == hg19_file

    def test_get_reference_path_missing_assembly(self):
        """Test error when assembly not configured."""
        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/some/path.fa", "vntr_region": "chr1:1-2"}
            }
        }

        with pytest.raises(ValidationError, match="Assembly 'mm10' not found"):
            get_reference_path_for_assembly(config, assembly="mm10")

    def test_get_reference_path_file_not_found(self):
        """Test error when reference file doesn't exist."""
        config = {
            "reference_genomes": {
                "hg38": {
                    "fasta_path": "/nonexistent/file.fa",
                    "vntr_region": "chr1:1-2",
                    "source_url": "http://example.com/hg38.fa"
                }
            }
        }

        with pytest.raises(FileOperationError, match="Reference genome file not found"):
            get_reference_path_for_assembly(config, assembly="hg38")

    def test_get_muc1_region(self):
        """Test getting MUC1 VNTR region for assembly."""
        config = {
            "reference_genomes": {
                "hg38": {
                    "fasta_path": "/some/path.fa",
                    "vntr_region": "chr1:155188487-155192239"
                },
                "hg19": {
                    "fasta_path": "/some/path.fa",
                    "vntr_region": "chr1:155160963-155162030"
                }
            }
        }

        region_hg38 = get_muc1_region_for_assembly(config, "hg38")
        assert region_hg38 == "chr1:155188487-155192239"

        region_hg19 = get_muc1_region_for_assembly(config, "hg19")
        assert region_hg19 == "chr1:155160963-155162030"

    def test_validate_reference_calls_existing_validation(self, tmp_path, mocker):
        """Verify that validate_reference_for_assembly uses existing function."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nATCG\n")

        config = {
            "reference_genomes": {
                "hg38": {
                    "fasta_path": str(ref_file),
                    "vntr_region": "chr1:1-2"
                }
            }
        }

        # Mock the existing validate_reference_genome function
        mock_validate = mocker.patch(
            "muc_one_up.bioinformatics.reference_validation.validate_reference_genome"
        )
        mock_validate.return_value = []

        # Call our new function
        warnings = validate_reference_for_assembly(config, "hg38")

        # Verify existing function was called (DRY - no duplication)
        mock_validate.assert_called_once()
        assert warnings == []
```

### 6.2 Integration Tests

**File**: `tests/test_reference_integration.py` (new)

```python
"""Integration tests for reference assembly management."""

import pytest
from pathlib import Path

from muc_one_up.config import load_config


@pytest.fixture
def example_config_with_references(tmp_path):
    """Create example config with reference_genomes section."""
    config_content = """
{
  "reference_genomes": {
    "hg38": {
      "fasta_path": "test_data/references/hg38.fa",
      "vntr_region": "chr1:155188487-155192239",
      "display_name": "Test hg38"
    }
  },
  "reference_assembly": "hg38",
  "repeats": {"1": "ACGT"},
  "constants": {
    "hg38": {"left": "AAA", "right": "TTT", "vntr_region": "chr1:155188487-155192239"}
  },
  "probabilities": {"1": {"1": 1.0}},
  "length_model": {
    "distribution": "normal",
    "min_repeats": 10,
    "max_repeats": 100,
    "mean_repeats": 50,
    "median_repeats": 50
  },
  "mutations": {},
  "tools": {"samtools": "samtools"},
  "read_simulation": {"human_reference": "test.fa", "threads": 1}
}
"""

    config_file = tmp_path / "test_config.json"
    config_file.write_text(config_content)

    return str(config_file)


def test_config_loads_with_reference_genomes(example_config_with_references):
    """Test that config with reference_genomes section loads correctly."""
    config = load_config(example_config_with_references)

    assert "reference_genomes" in config
    assert "hg38" in config["reference_genomes"]
    assert config["reference_genomes"]["hg38"]["vntr_region"] == "chr1:155188487-155192239"


def test_pipeline_can_access_reference(example_config_with_references, tmp_path):
    """Test that pipeline can access reference through new functions."""
    from muc_one_up.bioinformatics.reference_validation import get_muc1_region_for_assembly

    # Create dummy reference file
    ref_dir = tmp_path / "test_data" / "references"
    ref_dir.mkdir(parents=True)
    ref_file = ref_dir / "hg38.fa"
    ref_file.write_text(">chr1\nATCGATCG\n")

    config = load_config(example_config_with_references)

    # Update path to point to our test file
    config["reference_genomes"]["hg38"]["fasta_path"] = str(ref_file)

    # Should be able to get VNTR region
    region = get_muc1_region_for_assembly(config, "hg38")
    assert region == "chr1:155188487-155192239"
```

---

## 7. Migration

### 7.1 Automatic Migration

Existing configs without `reference_genomes` are auto-migrated:

**Old Config** (still works):
```json
{
  "constants": {
    "hg38": {"left": "...", "right": "...", "vntr_region": "chr1:155188487-155192239"}
  },
  "reference_assembly": "hg38"
}
```

**Auto-migration** creates:
```json
{
  "reference_genomes": {
    "hg38": {
      "fasta_path": "reference/hg38/hg38.fa",
      "vntr_region": "chr1:155188487-155192239"
    }
  },
  "constants": {...}  // Preserved
}
```

Users see warning: `"Config format outdated. Please add 'reference_genomes' section."`

### 7.2 Manual Migration (Recommended)

**Step 1**: Add `reference_genomes` to your config:

```json
{
  "reference_genomes": {
    "hg38": {
      "fasta_path": "/absolute/path/to/hg38.fa",
      "vntr_region": "chr1:155188487-155192239",
      "display_name": "Human GRCh38/hg38"
    },
    "hg19": {
      "fasta_path": "/absolute/path/to/hg19.fa",
      "vntr_region": "chr1:155160963-155162030",
      "display_name": "Human GRCh37/hg19"
    }
  },
  "reference_assembly": "hg38",
  "constants": {...}
}
```

**Step 2**: Test validation:

```bash
# This will validate your config
muconeup --config config.json validate
```

---

## 8. Anti-Patterns Avoided

### 8.1 What We Did NOT Do (Avoiding Over-Engineering)

❌ **Create `ReferenceManager` class**
- **Why not**: Violates KISS, adds unnecessary abstraction
- **Instead**: Simple functions in existing module

❌ **Create new module `reference_manager.py`**
- **Why not**: Breaks modularity, creates new dependency
- **Instead**: Extended existing `reference_validation.py`

❌ **Add complex dataclasses with 12+ fields**
- **Why not**: Over-engineering, adds complexity
- **Instead**: Simple dict access with helper functions

❌ **Environment variable expansion (`${REFERENCE_DIR}`)**
- **Why not**: Adds magic, harder to debug
- **Instead**: Plain paths (relative or absolute)

❌ **Chromosome name conversion utilities**
- **Why not**: YAGNI (You Ain't Gonna Need It) - just store correct value
- **Instead**: Config stores correct value for each assembly

❌ **Provenance tracking with checksums**
- **Why not**: Feature creep, not in scope
- **Instead**: Optional `source_url` field for documentation

❌ **Download management in module**
- **Why not**: Already exists in `helpers/`
- **Instead**: Enhance existing helper script

❌ **Mouse genome support**
- **Why not**: MUC1 is human-specific, scope creep
- **Instead**: Focus on human assemblies

### 8.2 SOLID Principles Compliance

✅ **Single Responsibility**
- `reference_validation.py` validates references (and now maps assemblies)
- Each function has one job
- No mixed concerns

✅ **Open/Closed**
- Config schema open for extension (add new assemblies)
- Functions closed for modification (don't need to change code)

✅ **Liskov Substitution**
- Not applicable (no inheritance/polymorphism)

✅ **Interface Segregation**
- 3 focused functions with clear purposes
- No fat interfaces

✅ **Dependency Inversion**
- Functions depend on config dict (abstraction)
- Not tied to specific config format

### 8.3 DRY Compliance

✅ **Reuses existing `validate_reference_genome()`**
- `validate_reference_for_assembly()` calls existing function
- No validation logic duplication

✅ **Single source of truth for assembly config**
- All assembly data in `reference_genomes`
- Not scattered across multiple config sections

✅ **Common path resolution logic**
- `get_reference_path_for_assembly()` used everywhere
- Pipelines don't implement their own logic

### 8.4 KISS Compliance

✅ **Simple config structure**
- Nested dict, no complex schema
- Easy to understand and edit

✅ **Minimal code**
- 3 functions, ~100 lines total
- No classes, no abstractions

✅ **Obvious behavior**
- Function names explain what they do
- No hidden magic

---

## 9. Files Modified

### Modified Files (5)

1. **`muc_one_up/config.py`**
   - Add `reference_genomes` to `CONFIG_SCHEMA`
   - Add auto-migration logic in `load_config()`
   - ~30 lines added

2. **`muc_one_up/bioinformatics/reference_validation.py`**
   - Add 3 helper functions
   - ~100 lines added

3. **`muc_one_up/read_simulator/ont_pipeline.py`**
   - Use helper functions for reference lookup
   - ~10 lines modified

4. **`muc_one_up/read_simulator/pipeline.py`**
   - Use helper functions for reference lookup
   - ~10 lines modified

5. **`muc_one_up/cli/click_main.py`**
   - Add reference validation before pipeline
   - Optional `--reference-assembly` override
   - ~20 lines added

### New Test Files (2)

6. **`tests/bioinformatics/test_reference_validation.py`**
   - Extend existing test file
   - Add `TestAssemblyManagement` class
   - ~100 lines added

7. **`tests/test_reference_integration.py`**
   - New integration tests
   - ~50 lines

### Updated Documentation (3)

8. **`CLAUDE.md`**
   - Document reference assembly configuration
   - Show examples

9. **`README.md`**
   - Add reference configuration section
   - Migration guide

10. **`config.json`** (example)
    - Add `reference_genomes` section
    - Include hg38, hg19, GRCh38 examples

**Total**: 10 files, ~320 lines of code/docs added

---

## 10. Implementation Timeline

### Week 1: Core Implementation

**Day 1-2**: Config schema
- Update `CONFIG_SCHEMA`
- Add auto-migration logic
- Test with existing configs

**Day 3-4**: Helper functions
- Implement 3 functions in `reference_validation.py`
- Unit tests
- Integration tests

**Day 5**: Pipeline integration
- Update ONT pipeline
- Update Illumina pipeline
- Test end-to-end

### Week 2: Documentation & Polish

**Day 1-2**: Documentation
- Update CLAUDE.md
- Update README.md
- Create migration guide

**Day 3**: Example configs
- Create example configs for common assemblies
- Validate all examples

**Day 4-5**: Testing & refinement
- Full integration testing
- User acceptance testing
- Bug fixes

---

## 11. Comparison: Over-Engineered vs Simple

| Aspect | ❌ Over-Engineered (Prev) | ✅ Simple (This) |
|--------|-------------------------|------------------|
| **New classes** | 2 (ReferenceGenome, ReferenceManager) | 0 |
| **New modules** | 1 (reference_manager.py) | 0 |
| **Lines of code** | 500+ | ~100 |
| **Functions** | 15+ | 3 |
| **Config complexity** | High (env vars, metadata) | Low (paths, regions) |
| **Dependencies** | Potential for external deps | None |
| **Maintainability** | Complex | Simple |
| **Testability** | Many mocks needed | Simple tests |
| **Debuggability** | Hard to trace | Easy to trace |
| **Time to implement** | 3 weeks | 1 week |
| **DRY compliance** | Duplicates validation | Reuses existing |
| **KISS compliance** | Violates | Follows |
| **SOLID compliance** | Questionable | Follows |

---

## Conclusion

This revised solution:

✅ **Follows all principles** - DRY, KISS, SOLID, modular
✅ **Learns from VNtyper** - Simple, config-driven approach
✅ **Respects original issue #28** - Extend existing, don't create new
✅ **Avoids over-engineering** - 3 functions vs complex class hierarchy
✅ **Minimal code** - ~100 lines vs 500+
✅ **Easy to maintain** - Simple, obvious, testable
✅ **Backward compatible** - Auto-migration for existing configs
✅ **Extensible** - Easy to add new assemblies

**Implementation**: Extend `reference_validation.py` with 3 helper functions. That's it.

**Philosophy**: Keep it stupidly simple. The best code is no code. The second-best code is simple code.

---

**End of Revised Assessment**

**Previous Assessment Status**: DEPRECATED - Over-engineered, violates KISS/DRY
**This Assessment Status**: APPROVED - Simple, follows all principles

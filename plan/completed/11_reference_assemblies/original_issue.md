# Issue #28: Human Reference Assembly Management (CORRECTED)

## CORRECTION NOTICE
**Original plan violated DRY principle by duplicating existing `reference_validation.py` module.**

Existing functionality (DO NOT DUPLICATE):
- `muc_one_up/bioinformatics/reference_validation.py`
  - `validate_reference_genome()` - checks FASTA + indices
  - `_validate_bwa_indices()` - checks BWA index files
  - `_validate_minimap2_indices()` - checks minimap2 index

## Problem
Config supports hg19/hg38 constants, but lacks helper functions to manage reference paths and validate assemblies. Users must manually track which reference file corresponds to which assembly.

## Current State
```json
// config.json
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
```

**Missing**: Link between assembly name and reference FASTA path

## Implementation (Extend Existing Module)

### 1. Config Schema Extension
```json
{
  "reference_genomes": {
    "hg38": {
      "fasta_path": "/path/to/hg38.fa",
      "assembly_name": "GRCh38"
    },
    "hg19": {
      "fasta_path": "/path/to/hg19.fa",
      "assembly_name": "GRCh37"
    }
  },
  "reference_assembly": "hg38"
}
```

### 2. Extend Existing reference_validation.py (MODIFY, DON'T CREATE NEW CLASS)
```python
# muc_one_up/bioinformatics/reference_validation.py
# ADD these functions to EXISTING module

from ..type_defs import ConfigDict

def get_reference_path_for_assembly(
    config: ConfigDict,
    assembly: str | None = None,
) -> Path:
    """Get reference FASTA path for specified assembly.

    Args:
        config: Configuration dictionary
        assembly: Assembly name (hg19/hg38). If None, uses config["reference_assembly"]

    Returns:
        Path to reference FASTA

    Raises:
        ValidationError: If assembly not configured
        FileOperationError: If reference file missing
    """
    # Get active assembly from config if not specified
    if assembly is None:
        assembly = config.get("reference_assembly", "hg38")

    # Get reference config for assembly
    ref_genomes = config.get("reference_genomes", {})
    ref_config = ref_genomes.get(assembly)

    if not ref_config:
        available = list(ref_genomes.keys())
        raise ValidationError(
            f"No reference configuration for assembly '{assembly}'. "
            f"Available: {available}"
        )

    # Get FASTA path
    fasta_path = ref_config.get("fasta_path")
    if not fasta_path:
        raise ValidationError(
            f"Missing 'fasta_path' for assembly '{assembly}' in config"
        )

    ref_path = Path(fasta_path).expanduser().resolve()

    # Validate using EXISTING function
    if not ref_path.exists():
        raise FileOperationError(f"Reference genome not found: {ref_path}")

    return ref_path


def validate_reference_for_assembly(
    config: ConfigDict,
    assembly: str | None = None,
    aligner: str = "bwa",
) -> list[str]:
    """Validate reference genome and indices for assembly.

    Uses EXISTING validate_reference_genome() function.

    Args:
        config: Configuration dictionary
        assembly: Assembly name (hg19/hg38)
        aligner: Aligner to check indices for

    Returns:
        List of warnings (empty if all OK)

    Raises:
        ValidationError: If assembly not configured
        FileOperationError: If reference missing
    """
    ref_path = get_reference_path_for_assembly(config, assembly)

    # Use EXISTING validation function - no duplication
    warnings = validate_reference_genome(ref_path, aligner=aligner)

    for warning in warnings:
        logging.warning(f"[{assembly}] {warning}")

    return warnings


def get_muc1_region_for_assembly(config: ConfigDict, assembly: str | None = None) -> str:
    """Get MUC1 VNTR region coordinates for assembly.

    Args:
        config: Configuration dictionary
        assembly: Assembly name

    Returns:
        Genomic region string (e.g., "chr1:155188487-155192239")
    """
    if assembly is None:
        assembly = config.get("reference_assembly", "hg38")

    constants = config.get("constants", {})
    assembly_constants = constants.get(assembly)

    if not assembly_constants:
        raise ValidationError(f"No constants defined for assembly: {assembly}")

    vntr_region = assembly_constants.get("vntr_region")
    if not vntr_region:
        raise ValidationError(f"No vntr_region defined for assembly: {assembly}")

    return vntr_region
```

### 3. Update Pipelines to Use Helper Functions (MODIFY)
```python
# muc_one_up/read_simulator/ont_pipeline.py

from ..bioinformatics.reference_validation import get_reference_path_for_assembly

def simulate_ont_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None
) -> str:
    # If human_reference not provided, get from config
    if human_reference is None:
        try:
            ref_path = get_reference_path_for_assembly(config)
            human_reference = str(ref_path)
            logging.info(f"Using reference: {human_reference}")
        except (ValidationError, FileOperationError) as e:
            logging.warning(f"Reference not configured: {e}")
            logging.warning("Will align to simulated reference instead")
            human_reference = input_fa

    # Rest of pipeline...
```

**Similarly update**: `pipeline.py` (Illumina)

### 4. Helper Script Enhancement (MODIFY)
```python
# helpers/download_references.py

from pathlib import Path
import urllib.request

ASSEMBLIES = {
    "hg38": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "checksum": "sha256:...",  # Add checksums
    },
    "hg19": {
        "url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        "checksum": "sha256:...",
    },
}

def download_reference(assembly: str, output_dir: Path) -> Path:
    """Download and prepare reference genome."""
    if assembly not in ASSEMBLIES:
        raise ValueError(f"Unknown assembly: {assembly}")

    config = ASSEMBLIES[assembly]
    output_path = output_dir / f"{assembly}.fa"

    # Download if missing
    if not output_path.exists():
        print(f"Downloading {assembly}...")
        urllib.request.urlretrieve(config["url"], str(output_path) + ".gz")
        # Decompress, verify checksum

    # Build indices
    build_bwa_index(output_path)
    build_samtools_index(output_path)

    return output_path
```

## Testing

### Unit Tests
```python
# tests/bioinformatics/test_reference_validation.py
# EXTEND existing test file

def test_get_reference_path_for_assembly(tmp_path):
    """Test reference path resolution."""
    ref_file = tmp_path / "hg38.fa"
    ref_file.write_text(">chr1\nATCG\n")

    config = {
        "reference_genomes": {
            "hg38": {"fasta_path": str(ref_file)}
        },
        "reference_assembly": "hg38"
    }

    path = get_reference_path_for_assembly(config)
    assert path == ref_file

def test_get_reference_path_missing_assembly():
    """Test error when assembly not configured."""
    config = {"reference_genomes": {}, "reference_assembly": "hg38"}

    with pytest.raises(ValidationError, match="No reference configuration"):
        get_reference_path_for_assembly(config)

def test_validate_reference_uses_existing_function(mocker, tmp_path):
    """Verify we use existing validate_reference_genome()."""
    ref_file = tmp_path / "hg38.fa"
    ref_file.write_text(">chr1\nATCG\n")

    config = {
        "reference_genomes": {"hg38": {"fasta_path": str(ref_file)}},
        "reference_assembly": "hg38"
    }

    # Mock the EXISTING function
    mock_validate = mocker.patch(
        "muc_one_up.bioinformatics.reference_validation.validate_reference_genome"
    )
    mock_validate.return_value = []

    validate_reference_for_assembly(config)

    # Verify existing function was called - no duplication
    mock_validate.assert_called_once()
```

## Documentation

### README.md
```markdown
### Reference Genome Configuration

Configure reference genomes in `config.json`:

```json
{
  "reference_genomes": {
    "hg38": {"fasta_path": "/data/hg38.fa"},
    "hg19": {"fasta_path": "/data/hg19.fa"}
  },
  "reference_assembly": "hg38"
}
```

**Download references:**
```bash
python helpers/download_references.py --assembly hg38 --output-dir references/
```

**Validate indices:**
```bash
muconeup --config config.json validate-reference
```
```

## Files Modified (NO NEW CLASSES)

### MODIFIED (3):
- `muc_one_up/bioinformatics/reference_validation.py` (add 3 helper functions)
- `muc_one_up/read_simulator/ont_pipeline.py` (use helper)
- `muc_one_up/read_simulator/pipeline.py` (use helper)
- `helpers/download_references.py` (enhance)
- `tests/bioinformatics/test_reference_validation.py` (extend tests)

### DO NOT CREATE:
- ❌ `muc_one_up/reference_manager.py` - Unnecessary abstraction
- ❌ `ReferenceManager` class - Violates KISS, duplicates validation

## Compliance with Programming Principles

✅ **DRY**: Extends existing module, doesn't duplicate validation logic

✅ **KISS**: Simple helper functions, no unnecessary classes

✅ **SOLID**: Single Responsibility maintained (validation module validates)

✅ **Modular**: Functions added to appropriate existing module

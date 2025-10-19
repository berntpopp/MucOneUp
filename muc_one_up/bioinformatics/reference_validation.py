"""Reference genome validation for MucOneUp.

This module validates reference genome files and their indices for
alignment tools like BWA and minimap2.

New in Issue #28:
- get_reference_path_for_assembly(): Get FASTA path for assembly
- validate_reference_for_assembly(): Validate reference for assembly
- get_muc1_region_for_assembly(): Get MUC1 VNTR region for assembly
"""

import logging
from pathlib import Path

from ..exceptions import FileOperationError, ValidationError
from ..type_defs import ConfigDict, FilePath


def validate_reference_genome(
    reference_path: FilePath,
    aligner: str = "bwa",
) -> list[str]:
    """Validate reference genome and indices exist.

    Args:
        reference_path: Path to reference FASTA
        aligner: Aligner to check indices for ('bwa' or 'minimap2')

    Returns:
        List of warnings (empty if all checks pass)

    Raises:
        FileOperationError: If reference not found
        ValidationError: If required indices missing or aligner unknown
    """
    ref_path = Path(reference_path)
    warnings: list[str] = []

    # Check reference exists
    if not ref_path.exists():
        raise FileOperationError(f"Reference genome not found: {reference_path}")

    # Check FASTA index (.fai)
    fai_path = Path(f"{reference_path}.fai")
    if not fai_path.exists():
        warnings.append(f"FASTA index missing: {fai_path}. Run: samtools faidx {reference_path}")

    # Check aligner-specific indices
    if aligner == "bwa":
        warnings.extend(_validate_bwa_indices(reference_path))
    elif aligner == "minimap2":
        warnings.extend(_validate_minimap2_indices(reference_path))
    else:
        raise ValidationError(f"Unknown aligner: '{aligner}'. Supported aligners: bwa, minimap2")

    return warnings


def _validate_bwa_indices(reference_path: FilePath) -> list[str]:
    """Validate BWA index files exist.

    Args:
        reference_path: Path to reference FASTA

    Returns:
        List of warnings for missing indices
    """
    required_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    missing = []

    for ext in required_extensions:
        index_path = Path(f"{reference_path}{ext}")
        if not index_path.exists():
            missing.append(ext)

    warnings: list[str] = []
    if missing:
        warnings.append(f"BWA index files missing: {missing}. Run: bwa index {reference_path}")

    return warnings


def _validate_minimap2_indices(reference_path: FilePath) -> list[str]:
    """Validate minimap2 index file exists.

    Args:
        reference_path: Path to reference FASTA

    Returns:
        List of warnings for missing indices
    """
    mmi_path = Path(f"{reference_path}.mmi")
    warnings: list[str] = []

    if not mmi_path.exists():
        warnings.append(
            f"Minimap2 index missing: {mmi_path}. Run: minimap2 -d {mmi_path} {reference_path}"
        )

    return warnings


def validate_bam_file(bam_path: FilePath, require_index: bool = True) -> list[str]:
    """Validate BAM file and optional index exist.

    Args:
        bam_path: Path to BAM file
        require_index: Require BAI index file (default: True)

    Returns:
        List of warnings

    Raises:
        FileOperationError: If BAM file not found
        ValidationError: If index required but missing
    """
    bam = Path(bam_path)
    warnings: list[str] = []

    # Check BAM exists
    if not bam.exists():
        raise FileOperationError(f"BAM file not found: {bam_path}")

    # Check for index
    bai_path = Path(f"{bam_path}.bai")
    if not bai_path.exists():
        msg = f"BAM index missing: {bai_path}. Run: samtools index {bam_path}"
        if require_index:
            raise ValidationError(msg)
        warnings.append(msg)

    return warnings


def validate_bed_file(bed_path: FilePath) -> None:
    """Validate BED file exists and has basic structure.

    Args:
        bed_path: Path to BED file

    Raises:
        FileOperationError: If file not found
        ValidationError: If BED format is invalid
    """
    bed = Path(bed_path)

    if not bed.exists():
        raise FileOperationError(f"BED file not found: {bed_path}")

    # Basic BED format validation
    with bed.open() as f:
        for i, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            # BED format: chrom start end [name] [score] [strand] ...
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValidationError(
                    f"Invalid BED format at line {i}: expected at least 3 fields, got {len(fields)}"
                )

            # Validate start/end are integers
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as e:
                raise ValidationError(f"Invalid BED coordinates at line {i}: {e}") from e

            # Validate start < end
            if start >= end:
                raise ValidationError(
                    f"Invalid BED coordinates at line {i}: start ({start}) >= end ({end})"
                )


# ============================================================================
# Issue #28: Reference Assembly Management
# ============================================================================


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
        raise ValidationError(f"No 'fasta_path' configured for assembly '{assembly}'")

    # Convert to Path and resolve
    fasta_path = Path(fasta_path_str)

    # If relative path, make it relative to current working directory
    if not fasta_path.is_absolute():
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


def get_muc1_region_for_assembly(config: ConfigDict, assembly: str | None = None) -> str:
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
        raise ValidationError("No 'reference_genomes' section in config")

    # Get assembly config
    assembly_config = ref_genomes.get(assembly)

    if not assembly_config:
        available = list(ref_genomes.keys())
        raise ValidationError(f"Assembly '{assembly}' not found. Available: {', '.join(available)}")

    # Get VNTR region
    vntr_region = assembly_config.get("vntr_region")

    if not vntr_region:
        raise ValidationError(f"No 'vntr_region' configured for assembly '{assembly}'")

    logging.debug(f"MUC1 VNTR region for {assembly}: {vntr_region}")
    return str(vntr_region)

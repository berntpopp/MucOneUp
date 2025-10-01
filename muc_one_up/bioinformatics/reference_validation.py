"""Reference genome validation for MucOneUp.

This module validates reference genome files and their indices for
alignment tools like BWA and minimap2.
"""

from pathlib import Path

from ..exceptions import FileOperationError, ValidationError
from ..type_defs import FilePath


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
        warnings.append(
            f"FASTA index missing: {fai_path}. " f"Run: samtools faidx {reference_path}"
        )

    # Check aligner-specific indices
    if aligner == "bwa":
        warnings.extend(_validate_bwa_indices(reference_path))
    elif aligner == "minimap2":
        warnings.extend(_validate_minimap2_indices(reference_path))
    else:
        raise ValidationError(
            f"Unknown aligner: '{aligner}'. " f"Supported aligners: bwa, minimap2"
        )

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
        warnings.append(f"BWA index files missing: {missing}. " f"Run: bwa index {reference_path}")

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
            f"Minimap2 index missing: {mmi_path}. " f"Run: minimap2 -d {mmi_path} {reference_path}"
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
                    f"Invalid BED format at line {i}: expected at least 3 fields, "
                    f"got {len(fields)}"
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

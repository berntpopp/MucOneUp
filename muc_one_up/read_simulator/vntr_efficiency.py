"""
VNTR Capture Efficiency Modeling.

Models probe capture bias in VNTR regions by selectively downsampling
reads to reflect real-world targeted sequencing behavior.

This module implements the empirically validated approach (October 2024)
using penalty factor 0.375 derived from comparison of 277 real Twist v2
samples vs 59 simulated samples.
"""

import logging
from pathlib import Path

from muc_one_up.read_simulator.utils import bed, samtools

logger = logging.getLogger(__name__)


# Default MUC1 VNTR region (hg38)
DEFAULT_VNTR_REGION = {
    "chr": "chr1",
    "start": 155188487,
    "end": 155192239,
    "name": "MUC1_VNTR",
}


class VNTREfficiencyError(Exception):
    """Exception raised for VNTR efficiency modeling failures."""

    pass


class VNTREfficiencyModel:
    """
    Model VNTR capture efficiency bias in BAM files.

    Applies empirically validated penalty factor to model the reduced
    capture efficiency of probe-based sequencing in VNTR regions.

    Attributes:
        penalty_factor: Fraction of VNTR reads to retain (default: 0.375)
        seed: Random seed for reproducible downsampling
        threads: Number of threads for samtools operations
        vntr_region: Dict defining VNTR coordinates
        capture_bed: Optional path to capture targets BED file
        flanking_size: Size of flanking regions if no capture BED (bp)
    """

    def __init__(
        self,
        penalty_factor: float = 0.375,
        seed: int = 42,
        threads: int = 8,
        vntr_region: dict | None = None,
        capture_bed: Path | None = None,
        flanking_size: int = 10000,
    ):
        """
        Initialize VNTR efficiency model.

        Args:
            penalty_factor: Fraction of VNTR reads to keep (0.1-1.0)
            seed: Random seed for reproducibility
            threads: Number of threads for samtools
            vntr_region: Dict with chr, start, end, name (default: MUC1 hg38)
            capture_bed: Path to capture targets BED (optional)
            flanking_size: Flanking region size if no capture BED

        Raises:
            ValueError: If parameters are invalid
            VNTREfficiencyError: If required tools not available
        """
        # Validate penalty factor
        if not 0.1 <= penalty_factor <= 1.0:
            raise ValueError(f"penalty_factor must be in [0.1, 1.0], got {penalty_factor}")

        self.penalty_factor = penalty_factor
        self.seed = seed
        self.threads = threads
        self.vntr_region = vntr_region or DEFAULT_VNTR_REGION
        self.capture_bed = Path(capture_bed) if capture_bed else None
        self.flanking_size = flanking_size

        # Validate required tools
        if not samtools.check_samtools_available():
            raise VNTREfficiencyError("samtools not found in PATH. Please install samtools.")

        logger.info("VNTREfficiencyModel initialized:")
        logger.info(f"  Penalty factor: {self.penalty_factor:.3f}")
        logger.info(f"  Seed: {self.seed}")
        logger.info(
            f"  VNTR region: {self.vntr_region['chr']}:"
            f"{self.vntr_region['start']}-{self.vntr_region['end']}"
        )

    def apply_efficiency_bias(self, input_bam: Path, output_bam: Path, temp_dir: Path) -> dict:
        """
        Apply VNTR capture efficiency bias to BAM file.

        Workflow:
        1. Create BED files for VNTR and non-VNTR regions
        2. Extract reads by region (mutually exclusive)
        3. Downsample VNTR reads by penalty factor
        4. Merge VNTR (downsampled) + non-VNTR (full)
        5. Sort and index output BAM
        6. Calculate coverage statistics

        Args:
            input_bam: Input BAM file (aligned reads)
            output_bam: Output BAM file (with efficiency bias)
            temp_dir: Directory for temporary files

        Returns:
            Dict with statistics:
                - vntr_coverage: Mean coverage in VNTR
                - non_vntr_coverage: Mean coverage in flanking/non-VNTR
                - coverage_ratio: VNTR/non-VNTR ratio
                - penalty_factor: Applied penalty factor
                - seed: Random seed used
                - input_reads: Total reads in input
                - output_reads: Total reads in output

        Raises:
            VNTREfficiencyError: If any step fails
            FileNotFoundError: If input_bam doesn't exist
        """
        # Validate inputs
        input_bam = Path(input_bam)
        output_bam = Path(output_bam)
        temp_dir = Path(temp_dir)

        if not input_bam.exists():
            raise FileNotFoundError(f"Input BAM not found: {input_bam}")

        if not samtools.validate_bam(input_bam):
            raise VNTREfficiencyError(f"Input BAM is corrupted: {input_bam}")

        # Create temp directory
        temp_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Applying VNTR efficiency bias to {input_bam.name}")
        logger.info(
            f"  Penalty factor: {self.penalty_factor:.3f} "
            f"(keep {self.penalty_factor * 100:.1f}% of VNTR reads)"
        )

        try:
            # Step 1: Create BED files
            vntr_bed, non_vntr_bed = self._create_bed_files(temp_dir)

            # Step 2: Extract reads by region
            vntr_reads_bam = temp_dir / "vntr_reads.bam"
            non_vntr_reads_bam = temp_dir / "non_vntr_reads.bam"

            logger.info("  Extracting VNTR reads...")
            self._extract_reads(input_bam, vntr_bed, vntr_reads_bam)

            logger.info("  Extracting non-VNTR reads...")
            self._extract_reads(input_bam, non_vntr_bed, non_vntr_reads_bam)

            # Step 3: Downsample VNTR reads
            vntr_downsampled_bam = temp_dir / "vntr_downsampled.bam"

            logger.info(f"  Downsampling VNTR reads by {self.penalty_factor:.3f}...")
            self._downsample_reads(vntr_reads_bam, vntr_downsampled_bam, self.penalty_factor)

            # Step 4: Merge and sort
            merged_bam = temp_dir / "merged.bam"

            logger.info("  Merging VNTR (downsampled) + non-VNTR (full)...")
            self._merge_and_sort(vntr_downsampled_bam, non_vntr_reads_bam, merged_bam, output_bam)

            # Step 5: Index output
            logger.info("  Indexing output BAM...")
            samtools.index_bam(output_bam)

            # Step 6: Calculate statistics
            logger.info("  Calculating coverage statistics...")
            stats = self._calculate_statistics(input_bam, output_bam, vntr_bed, non_vntr_bed)

            logger.info("  VNTR efficiency bias applied successfully")
            logger.info(f"    VNTR coverage: {stats['vntr_coverage']:.2f}x")
            logger.info(f"    Non-VNTR coverage: {stats['non_vntr_coverage']:.2f}x")
            logger.info(f"    Coverage ratio: {stats['coverage_ratio']:.3f}")

            return stats

        except Exception as e:
            logger.error(f"Failed to apply VNTR efficiency bias: {e}")
            raise VNTREfficiencyError(f"VNTR efficiency modeling failed: {e}") from e

    def _create_bed_files(self, temp_dir: Path) -> tuple[Path, Path]:
        """
        Create BED files for VNTR and non-VNTR regions.

        Args:
            temp_dir: Directory for BED files

        Returns:
            Tuple of (vntr_bed, non_vntr_bed)
        """
        logger.debug("Creating BED files...")

        # VNTR region BED
        vntr_bed = bed.create_vntr_bed(temp_dir, self.vntr_region)

        # Non-VNTR region BED
        if self.capture_bed and self.capture_bed.exists():
            # Use capture targets minus VNTR (preferred)
            logger.debug("Using capture BED for non-VNTR regions")
            non_vntr_bed = temp_dir / "non_vntr_region.bed"

            bed.create_non_vntr_bed_from_capture(non_vntr_bed, self.capture_bed, vntr_bed)
        else:
            # Fallback: Use flanking regions
            if self.capture_bed:
                logger.warning(
                    f"Capture BED not found: {self.capture_bed}, "
                    f"using flanking regions (±{self.flanking_size}bp)"
                )
            else:
                logger.debug(f"Using flanking regions (±{self.flanking_size}bp) for non-VNTR")

            _, _, non_vntr_bed = bed.create_flanking_beds(
                temp_dir, self.vntr_region, self.flanking_size
            )

        return vntr_bed, non_vntr_bed

    def _extract_reads(self, input_bam: Path, region_bed: Path, output_bam: Path) -> None:
        """
        Extract reads overlapping a region.

        Args:
            input_bam: Input BAM file
            region_bed: BED file defining region
            output_bam: Output BAM file
        """
        samtools.extract_reads_by_region(input_bam, region_bed, output_bam)

    def _downsample_reads(self, input_bam: Path, output_bam: Path, fraction: float) -> None:
        """
        Downsample reads by fraction.

        Args:
            input_bam: Input BAM file
            output_bam: Output BAM file
            fraction: Fraction of reads to keep
        """
        samtools.downsample_bam(input_bam, output_bam, fraction, self.seed)

    def _merge_and_sort(self, bam1: Path, bam2: Path, merged_bam: Path, sorted_bam: Path) -> None:
        """
        Merge two BAM files and sort by coordinate.

        Args:
            bam1: First BAM file (VNTR downsampled)
            bam2: Second BAM file (non-VNTR full)
            merged_bam: Output merged BAM (unsorted)
            sorted_bam: Output sorted BAM
        """
        # Merge
        samtools.merge_bams([bam1, bam2], merged_bam, threads=self.threads)

        # Sort
        samtools.sort_bam(merged_bam, sorted_bam, threads=self.threads)

    def _calculate_statistics(
        self, input_bam: Path, output_bam: Path, vntr_bed: Path, non_vntr_bed: Path
    ) -> dict:
        """
        Calculate coverage statistics for input and output BAMs.

        Args:
            input_bam: Input BAM file
            output_bam: Output BAM file
            vntr_bed: VNTR region BED
            non_vntr_bed: Non-VNTR region BED

        Returns:
            Dict with coverage statistics
        """
        # Read counts
        input_reads = samtools.get_bam_read_count(input_bam)
        output_reads = samtools.get_bam_read_count(output_bam)

        # Coverage in output BAM
        vntr_cov = samtools.calculate_mean_coverage(output_bam, vntr_bed)
        non_vntr_cov = samtools.calculate_mean_coverage(output_bam, non_vntr_bed)

        # Calculate ratio
        ratio = vntr_cov / non_vntr_cov if non_vntr_cov > 0 else 0.0

        return {
            "vntr_coverage": vntr_cov,
            "non_vntr_coverage": non_vntr_cov,
            "coverage_ratio": ratio,
            "penalty_factor": self.penalty_factor,
            "seed": self.seed,
            "input_reads": input_reads,
            "output_reads": output_reads,
            "reads_removed": input_reads - output_reads,
            "retention_fraction": output_reads / input_reads if input_reads > 0 else 0.0,
        }

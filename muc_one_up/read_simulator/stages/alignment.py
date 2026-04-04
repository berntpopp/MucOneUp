"""Illumina alignment and refinement stage (pipeline stages 9-11).

Extracts stages 9-11 of the Illumina read simulation pipeline into a
standalone function with a clean, testable interface.

Stages covered:
    9.  Align reads with bwa mem (align_reads wrapper)
    10. Apply VNTR capture efficiency bias (optional, enabled by default)
    11. Downsample to target coverage (optional)
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any

from ...exceptions import ConfigurationError
from ..assembly_context import AssemblyContext
from ..wrappers.bwa_wrapper import align_reads
from ..wrappers.samtools_wrapper import (
    calculate_target_coverage,
    calculate_vntr_coverage,
    downsample_bam,
    downsample_entire_bam,
)
from . import AlignmentResult

# Late import — resolved at call time to avoid circular imports:
#   from ..vntr_efficiency import VNTREfficiencyModel
# We import it at module level so that tests can patch it cleanly at
# ``muc_one_up.read_simulator.stages.alignment.VNTREfficiencyModel``.
try:
    from ..vntr_efficiency import VNTREfficiencyModel
except ImportError:  # pragma: no cover
    VNTREfficiencyModel = None  # type: ignore[assignment,misc]

logger = logging.getLogger(__name__)


def align_and_refine(
    tools: dict[str, str],
    rs_config: dict[str, Any],
    r1: str,
    r2: str,
    human_ref: str,
    output_dir: Path,
    output_base: str,
    assembly_ctx: AssemblyContext,
) -> AlignmentResult:
    """Run Illumina alignment and optional refinement (pipeline stages 9-11).

    Args:
        tools: Mapping of tool names to executable paths/commands.
        rs_config: Read-simulation configuration dictionary.
        r1: Path to the R1 paired-end FASTQ file.
        r2: Path to the R2 paired-end FASTQ file.
        human_ref: Path to the human reference FASTA.
        output_dir: Directory where output and intermediate files are written.
        output_base: Base name used for all output and intermediate file names.
        assembly_ctx: Resolved assembly context (provides vntr_region etc.).

    Returns:
        AlignmentResult with the final BAM path and lists of intermediate
        BAMs/files created during the stages.

    Raises:
        ConfigurationError: If downsampling is requested but required config
            values (``vntr_region`` or ``sample_target_bed``) are absent, or
            ``downsample_mode`` is invalid.
    """
    threads: int = rs_config.get("threads", 4)
    intermediate_bams: list[str] = []
    intermediate_files: list[str] = []

    # ------------------------------------------------------------------
    # Stage 9: Align reads with bwa mem
    # ------------------------------------------------------------------
    aligned_bam = str(output_dir / f"{output_base}.bam")
    logger.info("9. Aligning reads to reference: %s", human_ref)
    align_reads(r1, r2, human_ref, aligned_bam, tools, threads)
    current_bam = aligned_bam

    # ------------------------------------------------------------------
    # Stage 10: Apply VNTR capture efficiency bias (optional)
    # ------------------------------------------------------------------
    vntr_config = rs_config.get("vntr_capture_efficiency", {})
    vntr_enabled = vntr_config.get("enabled", True)  # enabled by default

    if vntr_enabled:
        logger.info("10. Applying VNTR capture efficiency bias")
        try:
            penalty_factor: float = vntr_config.get("penalty_factor", 0.375)
            vntr_seed: int = vntr_config.get("seed", rs_config.get("seed", 42))
            vntr_region = vntr_config.get("vntr_region")
            capture_bed = vntr_config.get("capture_bed")
            flanking_size: int = vntr_config.get("flanking_size", 10000)

            vntr_model = VNTREfficiencyModel(
                penalty_factor=penalty_factor,
                seed=vntr_seed,
                threads=threads,
                vntr_region=vntr_region,
                capture_bed=Path(capture_bed) if capture_bed else None,
                flanking_size=flanking_size,
            )

            # Track pre-bias BAM as intermediate before replacement
            intermediate_bams.append(current_bam)

            vntr_biased_bam = str(output_dir / f"{output_base}_vntr_biased.bam")
            temp_dir = output_dir / "_vntr_efficiency_temp"

            vntr_stats = vntr_model.apply_efficiency_bias(
                input_bam=Path(current_bam),
                output_bam=Path(vntr_biased_bam),
                temp_dir=temp_dir,
            )

            # Optionally save statistics JSON
            if vntr_config.get("validation", {}).get("report_statistics", True):
                stats_file = str(output_dir / f"{output_base}_vntr_efficiency_stats.json")
                with open(stats_file, "w") as fh:
                    json.dump(vntr_stats, fh, indent=2)
                logger.info("  VNTR efficiency statistics saved to %s", stats_file)

            current_bam = vntr_biased_bam
            logger.info("  VNTR efficiency bias applied successfully")

            # Optionally generate FASTQ from the VNTR-biased BAM
            vntr_fastq_config = vntr_config.get("output_fastq", {})
            vntr_fastq_enabled: bool = vntr_fastq_config.get("enabled", True)

            if vntr_fastq_enabled:
                logger.info("  Generating FASTQ files from VNTR-biased BAM")
                try:
                    from ..wrappers.samtools_wrapper import (
                        FastqConversionOptions,
                        convert_bam_to_paired_fastq,
                    )

                    vntr_base = Path(vntr_biased_bam).stem
                    vntr_fq1 = str(output_dir / f"{vntr_base}_R1.fastq.gz")
                    vntr_fq2 = str(output_dir / f"{vntr_base}_R2.fastq.gz")

                    opts = FastqConversionOptions(
                        output_singleton=vntr_fastq_config.get("singleton_file"),
                        preserve_read_names=vntr_fastq_config.get("preserve_read_names", True),
                        validate_pairs=True,
                        threads=threads,
                        timeout=1800,
                    )

                    fq1_out, fq2_out = convert_bam_to_paired_fastq(
                        samtools_cmd=tools["samtools"],
                        input_bam=vntr_biased_bam,
                        output_fq1=vntr_fq1,
                        output_fq2=vntr_fq2,
                        options=opts,
                    )
                    logger.info("  VNTR-biased FASTQ R1: %s", Path(fq1_out).name)
                    logger.info("  VNTR-biased FASTQ R2: %s", Path(fq2_out).name)
                except Exception as fastq_err:
                    logger.warning("  Failed to generate VNTR-biased FASTQ files: %s", fastq_err)
                    logger.warning("  Continuing with BAM output only")
            else:
                logger.info("  Skipping VNTR-biased FASTQ generation (disabled in config)")

            # Clean up temporary directory
            if temp_dir.exists():
                shutil.rmtree(temp_dir)

        except Exception as exc:
            logger.error("VNTR efficiency modeling failed: %s", exc)
            logger.warning("Continuing with original BAM file (no bias applied)")
            # Undo the intermediate_bams append if it was done before the failure
            if intermediate_bams and intermediate_bams[-1] == aligned_bam:
                intermediate_bams.pop()
            current_bam = aligned_bam
    else:
        logger.info("10. Skipping VNTR capture efficiency (disabled in config)")

    # ------------------------------------------------------------------
    # Stage 11: Optionally downsample to target coverage
    # ------------------------------------------------------------------
    target_coverage = rs_config.get("coverage")
    if target_coverage:
        logger.info("11. Downsampling to target coverage")
        mode: str = rs_config.get("downsample_mode", "vntr").strip().lower()

        if mode == "vntr":
            vntr_region = assembly_ctx.vntr_region
            if not vntr_region:
                raise ConfigurationError(
                    f"VNTR region not specified in config for {assembly_ctx.assembly_name}. "
                    f"Add 'vntr_region_{assembly_ctx.assembly_name}' to config"
                )
            current_cov, depth_file = calculate_vntr_coverage(
                tools["samtools"],
                current_bam,
                vntr_region,
                threads,
                str(output_dir),
                output_base,
            )
            region_info = vntr_region
        elif mode == "non_vntr":
            bed_file = rs_config.get("sample_target_bed")
            if not bed_file:
                raise ConfigurationError(
                    "For non-VNTR downsampling, 'sample_target_bed' must be provided in config"
                )
            current_cov, depth_file = calculate_target_coverage(
                tools["samtools"],
                current_bam,
                bed_file,
                threads,
                str(output_dir),
                output_base,
            )
            region_info = f"BED file: {bed_file}"
        else:
            raise ConfigurationError(
                f"Invalid downsample_mode '{mode}' in config; use 'vntr' or 'non_vntr'"
            )

        intermediate_files.append(depth_file)

        if current_cov > target_coverage:
            fraction = min(max(target_coverage / current_cov, 0.0), 1.0)
            logger.info(
                "Downsampling BAM from %.2fx to %.2fx (fraction: %.4f) based on %s",
                current_cov,
                target_coverage,
                fraction,
                region_info,
            )
            intermediate_bams.append(current_bam)
            downsampled_bam = str(output_dir / f"{output_base}_downsampled.bam")

            if mode == "vntr":
                downsample_bam(
                    tools["samtools"],
                    current_bam,
                    downsampled_bam,
                    vntr_region,
                    fraction,
                    rs_config.get("downsample_seed", 42),
                    threads,
                )
            else:
                downsample_entire_bam(
                    tools["samtools"],
                    current_bam,
                    downsampled_bam,
                    fraction,
                    rs_config.get("downsample_seed", 42),
                    threads,
                )
            current_bam = downsampled_bam
        else:
            logger.info(
                "Current coverage (%.2fx) is below the target; no downsampling performed.",
                current_cov,
            )
    else:
        logger.info("11. Skipping downsampling (no target coverage specified)")

    return AlignmentResult(
        final_bam=current_bam,
        intermediate_bams=intermediate_bams,
        intermediate_files=intermediate_files,
    )

"""Read simulation pipeline orchestrator.

Thin orchestrator that delegates to stage modules:
- stages.fragment_preparation: Stages 1-8
- stages.alignment: Stages 9-11
- stages.source_manifest: Stage 12a
- pipeline_utils: Shared utilities
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .output_config import OutputConfig

from .assembly_context import AssemblyContext
from .pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    resolve_human_reference,
    resolve_pipeline_outputs,
)
from .stages.alignment import align_and_refine
from .stages.fragment_preparation import prepare_fragments
from .stages.source_manifest import generate_read_manifest
from .utils import capture_tool_versions, check_external_tools, log_tool_versions

logger = logging.getLogger(__name__)


def simulate_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Run the complete Illumina read simulation pipeline.

    Orchestrates the multi-step process of simulating Illumina reads from a
    FASTA file by delegating to extracted stage modules:

    - Stages 1-8  (fragment preparation): replace Ns, systematic errors,
      2bit conversion, subset reference extraction, pblat alignment,
      fragment simulation, read creation, and FASTQ splitting.
    - Stages 9-11 (alignment and refinement): BWA MEM alignment, optional
      VNTR capture efficiency bias, and optional coverage downsampling.
    - Stage 12a   (source manifest): read provenance manifest (when
      source_tracker is provided).
    - Stage 12b   (metadata): pipeline metadata TSV.
    - Stage 12c   (cleanup): removal of intermediate files (configurable).

    Args:
        config: Dictionary containing "tools" and "read_simulation" sections.
        input_fa: Input simulated FASTA file.
        source_tracker: Optional read source tracker for provenance.
        output_config: Optional OutputConfig controlling output directory and
            base name. When provided, takes precedence over config-derived
            paths and input-file-derived naming.

    Returns:
        Path to the final output BAM file.

    Raises:
        ConfigurationError: If required configuration is missing.
    """
    start_time = datetime.now()
    logger.info(
        "Starting read simulation pipeline at %s",
        start_time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})
    assembly_ctx = AssemblyContext.from_configs(config, rs_config)

    # Validate and log tools
    check_external_tools(tools)
    logger.info("=" * 60)
    logger.info("Tool Version Information")
    logger.info("=" * 60)
    tools_to_check = ["reseq", "faToTwoBit", "samtools", "pblat", "bwa"]
    tools_subset = {k: v for k, v in tools.items() if k in tools_to_check}
    tool_versions = capture_tool_versions(tools_subset)
    log_tool_versions(tool_versions)
    logger.info("=" * 60)

    # Resolve output paths
    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs_config, output_config)

    # Legacy FASTQ path overrides (only when no output_config)
    if output_config is not None:
        reads_fq1 = str(output_dir / f"{output_base}_R1.fastq.gz")
        reads_fq2 = str(output_dir / f"{output_base}_R2.fastq.gz")
    else:
        reads_fq1 = rs_config.get("output_fastq1", str(output_dir / f"{output_base}_R1.fastq.gz"))
        reads_fq2 = rs_config.get("output_fastq2", str(output_dir / f"{output_base}_R2.fastq.gz"))

    logger.info("Output filenames:")
    logger.info("  FASTQ pair: %s, %s", reads_fq1, reads_fq2)

    # Stages 1-8: Fragment preparation
    frag_result = prepare_fragments(
        tools=tools,
        rs_config=rs_config,
        input_fa=input_fa,
        output_dir=output_dir,
        output_base=output_base,
        assembly_ctx=assembly_ctx,
        source_tracker=source_tracker,
    )

    # Resolve human reference
    human_ref = resolve_human_reference(config, assembly_ctx, aligner="bwa")

    # Stages 9-11: Alignment and refinement
    align_result = align_and_refine(
        tools=tools,
        rs_config=rs_config,
        r1=frag_result.r1_fastq,
        r2=frag_result.r2_fastq,
        human_ref=human_ref,
        output_dir=output_dir,
        output_base=output_base,
        assembly_ctx=assembly_ctx,
    )

    end_time = datetime.now()
    duration = end_time - start_time
    logger.info(
        "Pipeline completed at %s (duration: %s)",
        end_time.strftime("%Y-%m-%d %H:%M:%S"),
        str(duration).split(".")[0],
    )
    logger.info("Final outputs:")
    logger.info("  Aligned and indexed BAM: %s", align_result.final_bam)
    logger.info("  Paired FASTQ files: %s and %s", reads_fq1, reads_fq2)

    # Stage 12a: Source tracking manifest
    fragment_origins_path = (
        str(output_dir / f"{output_base}_fragment_origins.tsv")
        if source_tracker is not None
        else None
    )
    generate_read_manifest(
        sidecar_path=fragment_origins_path,
        final_bam=align_result.final_bam,
        input_fa=input_fa,
        source_tracker=source_tracker,
        output_dir=output_dir,
        output_base=output_base,
    )

    # Stage 12b: Metadata
    create_pipeline_metadata(
        output_dir=output_dir,
        output_base=output_base,
        config=config,
        start_time=start_time,
        end_time=end_time,
        platform="Illumina",
        tools_used=["reseq", "faToTwoBit", "pblat", "bwa", "samtools"],
    )

    # Stage 12c: Cleanup
    keep_intermediates = rs_config.get("keep_intermediate_files", False)
    if not keep_intermediates:
        all_intermediates: list[str | None] = list(frag_result.intermediate_files)
        for bam in align_result.intermediate_bams:
            all_intermediates.append(bam)
            bam_index = f"{bam}.bai"
            if Path(bam_index).exists():
                all_intermediates.append(bam_index)
        all_intermediates.extend(align_result.intermediate_files)

        final_bam_index = f"{align_result.final_bam}.bai"
        all_intermediates = [
            f for f in all_intermediates if f not in (align_result.final_bam, final_bam_index)
        ]

        logger.info("=" * 60)
        logger.info("Intermediate File Cleanup")
        logger.info("=" * 60)
        logger.info("Total files to remove: %d", len(all_intermediates))
        logger.info("Preserving final output: %s", Path(align_result.final_bam).name)
        cleanup_intermediates(all_intermediates)
        logger.info("Cleanup completed successfully")
        logger.info("=" * 60)
    else:
        logger.info("=" * 60)
        logger.info("Keeping intermediate files (keep_intermediate_files=true)")
        logger.info("=" * 60)

    return align_result.final_bam


# Backward compatibility CLI
if __name__ == "__main__":  # OK: top-level entry point
    import json
    import sys

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")

    if len(sys.argv) != 3:
        print("Usage: python pipeline.py <config.json> <input_fasta>")
        sys.exit(1)

    with Path(sys.argv[1]).open() as fh:
        cfg = json.load(fh)

    simulate_reads_pipeline(cfg, sys.argv[2])

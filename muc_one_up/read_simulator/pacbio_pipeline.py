#!/usr/bin/env python3
"""
PacBio HiFi read simulation pipeline.

This module orchestrates the complete PacBio HiFi read simulation workflow,
integrating pbsim3 for CLR simulation, CCS for consensus generation, and
minimap2 for alignment.

Key Functions:
    simulate_pacbio_hifi_reads: Complete HiFi simulation pipeline

Pipeline Stages:
    1. Multi-pass CLR simulation (pbsim3)
    2. HiFi consensus generation (CCS)
    3. BAM→FASTQ conversion (samtools)
    4. Read alignment (minimap2 with map-hifi preset)
    5. Intermediate file cleanup

Example:
    Simulate PacBio HiFi reads::

        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        final_bam = simulate_pacbio_hifi_reads(
            config={
                "tools": {
                    "pbsim3": "pbsim",
                    "ccs": "ccs",
                    "samtools": "samtools",
                    "minimap2": "minimap2"
                },
                "pacbio_params": {
                    "model_type": "errhmm",
                    "model_file": "/models/ERRHMM-SEQUEL.model",
                    "coverage": 30,
                    "pass_num": 3,
                    "min_passes": 3,
                    "min_rq": 0.99,
                    "threads": 8
                }
            },
            input_fa="vntr_diploid.fa",
            output_dir="/output",
            output_base="sample",
            human_reference="/ref/hg38.fa"
        )

Notes:
    - Automatically manages intermediate files
    - Validates all inputs before processing
    - Provides detailed progress logging
    - Handles errors gracefully with cleanup

References:
    - pbsim3: https://github.com/yukiteruono/pbsim3
    - CCS: https://ccs.how/
    - minimap2: https://github.com/lh3/minimap2
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .output_config import OutputConfig

from .constants import MINIMAP2_PRESET_PACBIO_HIFI
from .pipeline_utils import (
    cleanup_unless_kept,
    create_pipeline_metadata,
    resolve_pipeline_outputs,
)
from .utils.reference_utils import get_reference_info
from .wrappers.ccs_wrapper import run_ccs_consensus, validate_ccs_parameters
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_simulation, validate_pbsim3_parameters
from .wrappers.samtools_wrapper import convert_bam_to_fastq, merge_bam_files

__all__ = ["simulate_pacbio_hifi_reads"]


def _register_ccs_intermediates(bam_path: str, intermediate_files: list[str | None]) -> None:
    """Register CCS sidecar files alongside a HiFi BAM for cleanup."""
    base = bam_path.removesuffix(".bam")
    intermediate_files.extend(
        [
            f"{bam_path}.pbi",
            f"{bam_path}.gzi",
            f"{base}.ccs_report.txt",
            f"{base}.zmw_metrics.json.gz",
        ]
    )


def _simulate_haplotype_hifi(
    clr_bam: str,
    ccs_cmd: str,
    output_path: str,
    min_passes: int,
    min_rq: float,
    threads: int,
    seed: int | None,
    hap_idx: int,
) -> str:
    """Run CCS consensus on a single haplotype's CLR BAM.

    Args:
        clr_bam: Path to CLR BAM for this haplotype.
        ccs_cmd: CCS command path.
        output_path: Output HiFi BAM path.
        min_passes: Minimum CCS passes.
        min_rq: Minimum read quality.
        threads: Thread count.
        seed: Base seed (hap_idx added for uniqueness).
        hap_idx: 1-based haplotype index.

    Returns:
        Path to HiFi BAM.
    """
    hifi_seed = seed + hap_idx if seed is not None else None
    hifi_bam = run_ccs_consensus(
        ccs_cmd=ccs_cmd,
        input_bam=clr_bam,
        output_bam=output_path,
        min_passes=min_passes,
        min_rq=min_rq,
        threads=threads,
        seed=hifi_seed,
    )
    logging.info("Haplotype %d CCS complete: %s", hap_idx, hifi_bam)
    return hifi_bam


def simulate_pacbio_hifi_reads(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """
    Complete PacBio HiFi read simulation pipeline.

    Orchestrates multi-pass CLR simulation, HiFi consensus generation,
    FASTQ conversion, and alignment to produce publication-ready HiFi reads.
    """
    from typing import cast

    from ..type_defs import PacbioConfig

    pacbio_params = cast(PacbioConfig, config.get("pacbio_params", {}))
    tools = config.get("tools", {})

    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    ccs_cmd = tools.get("ccs", "ccs")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    model_type = pacbio_params.get("model_type")
    model_file = pacbio_params.get("model_file")
    coverage = pacbio_params.get("coverage")
    pass_num = pacbio_params.get("pass_num")
    min_passes = pacbio_params.get("min_passes")
    min_rq = pacbio_params.get("min_rq")
    threads = pacbio_params.get("threads", 4)
    seed = pacbio_params.get("seed")

    from .constants import (
        DEFAULT_PBSIM3_ACCURACY_MEAN,
        DEFAULT_PBSIM3_ACCURACY_MIN,
        DEFAULT_PBSIM3_ACCURACY_SD,
        DEFAULT_PBSIM3_LENGTH_MAX,
        DEFAULT_PBSIM3_LENGTH_MEAN,
        DEFAULT_PBSIM3_LENGTH_MIN,
        DEFAULT_PBSIM3_LENGTH_SD,
    )

    accuracy_mean = pacbio_params.get("accuracy_mean", DEFAULT_PBSIM3_ACCURACY_MEAN)
    accuracy_sd = pacbio_params.get("accuracy_sd", DEFAULT_PBSIM3_ACCURACY_SD)
    accuracy_min = pacbio_params.get("accuracy_min", DEFAULT_PBSIM3_ACCURACY_MIN)
    length_mean = pacbio_params.get("length_mean", DEFAULT_PBSIM3_LENGTH_MEAN)
    length_sd = pacbio_params.get("length_sd", DEFAULT_PBSIM3_LENGTH_SD)
    length_min = pacbio_params.get("length_min", DEFAULT_PBSIM3_LENGTH_MIN)
    length_max = pacbio_params.get("length_max", DEFAULT_PBSIM3_LENGTH_MAX)

    required_params = {
        "model_type": model_type,
        "model_file": model_file,
        "coverage": coverage,
        "pass_num": pass_num,
        "min_passes": min_passes,
        "min_rq": min_rq,
    }

    missing_params = [key for key, value in required_params.items() if value is None]
    if missing_params:
        raise RuntimeError(
            f"Missing required PacBio parameters in config: {', '.join(missing_params)}"
        )

    if model_type is None or model_file is None or coverage is None:
        raise RuntimeError("Required PacBio parameters missing")  # unreachable
    if pass_num is None or min_passes is None or min_rq is None:
        raise RuntimeError("Required PacBio parameters missing")  # unreachable

    try:
        validate_pbsim3_parameters(
            coverage=coverage,
            pass_num=pass_num,
            accuracy_mean=accuracy_mean,
            accuracy_sd=accuracy_sd,
            accuracy_min=accuracy_min,
            length_mean=length_mean,
            length_min=length_min,
            length_max=length_max,
        )

        validate_ccs_parameters(
            min_passes=min_passes,
            min_rq=min_rq,
            threads=threads,
        )
    except ValueError as e:
        raise RuntimeError(f"Invalid PacBio parameters: {e}") from e

    output_dir_path, output_base = resolve_pipeline_outputs(
        input_fa, config.get("read_simulation", {}), output_config
    )

    clr_prefix = str(output_dir_path / f"{output_base}_clr")
    clr_bam = f"{clr_prefix}.bam"

    hifi_bam = str(output_dir_path / f"{output_base}_hifi.bam")
    hifi_fastq = str(output_dir_path / f"{output_base}_hifi.fastq")
    aligned_bam = str(output_dir_path / f"{output_base}_aligned.bam")

    intermediate_files: list[str | None] = []
    start_time = datetime.now()

    try:
        logging.info("=" * 80)
        logging.info("STAGE 1: Multi-pass CLR simulation with pbsim3")
        logging.info("=" * 80)

        ref_info = get_reference_info(input_fa)

        # Register expected Stage 1 paths before pbsim3 runs so partial outputs are cleaned on failure
        expected_bases = [f"{clr_prefix}_{i:04d}" for i in range(1, ref_info.num_sequences + 1)]
        if ref_info.num_sequences == 1:
            expected_bases.append(clr_prefix)
        for base in expected_bases:
            intermediate_files.extend(
                [f"{base}.bam", f"{base}.sam", f"{base}.maf", f"{base}.maf.gz", f"{base}.ref"]
            )

        clr_bams = run_pbsim3_simulation(
            pbsim3_cmd=pbsim3_cmd,
            samtools_cmd=samtools_cmd,
            reference=input_fa,
            model_type=model_type,
            model_file=model_file,
            coverage=coverage / ref_info.num_sequences,
            output_prefix=clr_prefix,
            pass_num=pass_num,
            accuracy_mean=accuracy_mean,
            accuracy_sd=accuracy_sd,
            accuracy_min=accuracy_min,
            length_mean=length_mean,
            length_sd=length_sd,
            length_min=length_min,
            length_max=length_max,
            seed=seed,
            threads=threads,
        )

        intermediate_files.extend(clr_bams)
        maf_files_by_haplotype: dict[int, list[str]] = {}
        for hap_idx, clr_bam in enumerate(clr_bams, 1):
            base = clr_bam.removesuffix(".bam")
            maf_path = f"{base}.maf"
            maf_gz_path = f"{base}.maf.gz"
            intermediate_files.extend([maf_path, maf_gz_path, f"{base}.ref"])
            maf_files_by_haplotype[hap_idx] = [maf_path, maf_gz_path]

        logging.info("=" * 80)
        logging.info("STAGE 2: HiFi consensus generation with CCS")
        logging.info("=" * 80)

        intermediate_files.append(hifi_bam)
        _register_ccs_intermediates(hifi_bam, intermediate_files)

        hifi_bams = []
        for i, clr_bam in enumerate(clr_bams, 1):
            if len(clr_bams) > 1:
                hifi_bam_individual = str(output_dir_path / f"{output_base}_hifi_{i:04d}.bam")
                intermediate_files.append(hifi_bam_individual)
                _register_ccs_intermediates(hifi_bam_individual, intermediate_files)
                logging.info("Processing haplotype %d/%d: %s", i, len(clr_bams), clr_bam)
            else:
                hifi_bam_individual = hifi_bam

            hifi_bam_result = _simulate_haplotype_hifi(
                clr_bam=clr_bam,
                ccs_cmd=ccs_cmd,
                output_path=hifi_bam_individual,
                min_passes=min_passes,
                min_rq=min_rq,
                threads=threads,
                seed=seed,
                hap_idx=i,
            )
            hifi_bams.append(hifi_bam_result)

        if len(hifi_bams) > 1:
            logging.info(f"Merging {len(hifi_bams)} HiFi BAMs into: {hifi_bam}")
            hifi_bam = merge_bam_files(
                samtools_cmd=samtools_cmd,
                input_bams=hifi_bams,
                output_bam=hifi_bam,
                threads=threads,
            )
        else:
            hifi_bam = hifi_bams[0]

        logging.info("=" * 80)
        logging.info("STAGE 3: Converting HiFi BAM to FASTQ")
        logging.info("=" * 80)

        hifi_fastq = convert_bam_to_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=hifi_bam,
            output_fastq=hifi_fastq,
            threads=threads,
        )

        if source_tracker is not None:
            from .parsers.pacbio_parser import parse_pacbio_reads, reconcile_ccs_reads

            logging.info("Generating read source tracking manifest...")
            all_origins = []
            for hap_idx, maf_paths in maf_files_by_haplotype.items():
                existing_mafs = [p for p in maf_paths if Path(p).exists()]
                if existing_mafs:
                    origins = parse_pacbio_reads(
                        maf_paths=existing_mafs,
                        haplotype_index=hap_idx,
                    )
                    all_origins.extend(origins)

            all_origins = reconcile_ccs_reads(all_origins, hifi_fastq)

            manifest_path = f"{output_dir_path / output_base!s}_read_manifest.tsv.gz"
            source_tracker.write_manifest(source_tracker.annotate_reads(all_origins), manifest_path)
            logging.info(
                "Read source manifest written: %s (%d reads)", manifest_path, len(all_origins)
            )

        if human_reference is not None:
            intermediate_files.append(hifi_fastq)

            logging.info("=" * 80)
            logging.info("STAGE 4: Aligning HiFi reads with minimap2 (map-hifi preset)")
            logging.info("=" * 80)

            aligned_bam = align_reads_with_minimap2(
                minimap2_cmd=minimap2_cmd,
                samtools_cmd=samtools_cmd,
                reference=human_reference,
                reads_fastq=hifi_fastq,
                output_bam=aligned_bam,
                preset=MINIMAP2_PRESET_PACBIO_HIFI,
                threads=threads,
            )
            final_output = aligned_bam
        else:
            final_output = hifi_fastq

        end_time = datetime.now()
        duration = end_time - start_time

        logging.info("=" * 80)
        logging.info(
            "PacBio HiFi simulation pipeline complete! (duration: %s)",
            str(duration).split(".")[0],
        )
        logging.info(f"Final output: {final_output}")

        tools_used = ["pbsim3", "ccs", "samtools"]
        if human_reference is not None:
            tools_used.append("minimap2")

        metadata_file = create_pipeline_metadata(
            output_dir=output_dir_path,
            output_base=output_base,
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="PacBio",
            tools_used=tools_used,
        )
        logging.info("Metadata file: %s", metadata_file)
        logging.info("=" * 80)

        return final_output

    finally:
        try:
            if "intermediate_files" in locals():
                logging.info("=" * 80)
                logging.info("STAGE 5: Cleaning up intermediate files")
                logging.info("=" * 80)
                cleanup_unless_kept(config, intermediate_files)
        except Exception as cleanup_error:
            logging.warning("Cleanup failed (non-fatal): %s", cleanup_error)

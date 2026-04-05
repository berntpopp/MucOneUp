"""Amplicon read simulation pipeline for PacBio.

Orchestrates the complete amplicon simulation workflow:
1. Haplotype extraction (diploid) or pass-through (haploid)
2. Primer-based amplicon extraction per haplotype
3. PCR bias coverage split computation
4. Template FASTA generation (N copies per allele)
5. PBSIM3 template mode simulation per allele
6. CCS consensus generation
7. BAM merge (diploid)
8. Alignment to human reference (optional)
"""

from __future__ import annotations

import logging
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .output_config import OutputConfig

from ..exceptions import ExternalToolError, FileOperationError, ReadSimulationError
from .amplicon_common import extract_and_prepare_amplicons
from .constants import MINIMAP2_PRESET_PACBIO_HIFI
from .pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    resolve_pipeline_outputs,
)
from .wrappers.ccs_wrapper import run_ccs_consensus
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_wrapper import convert_bam_to_fastq, merge_bam_files


def _simulate_haplotype_amplicon(
    clr_bam: str,
    ccs_cmd: str,
    output_path: str,
    min_passes: int,
    min_rq: float,
    threads: int,
    seed: int | None,
    hap_idx: int,
    bam_idx: int,
) -> str:
    """Run CCS consensus on a single amplicon CLR BAM.

    Args:
        clr_bam: Path to CLR BAM.
        ccs_cmd: CCS command path.
        output_path: Output HiFi BAM path.
        min_passes: Minimum CCS passes.
        min_rq: Minimum read quality.
        threads: Thread count.
        seed: Base seed (hap_idx*100 + bam_idx added for uniqueness).
        hap_idx: 1-based haplotype index.
        bam_idx: 1-based BAM index within haplotype.

    Returns:
        Path to HiFi BAM.
    """
    hifi_seed = (seed + hap_idx * 100 + bam_idx) if seed is not None else None
    return run_ccs_consensus(
        ccs_cmd=ccs_cmd,
        input_bam=clr_bam,
        output_bam=output_path,
        min_passes=min_passes,
        min_rq=min_rq,
        threads=threads,
        seed=hifi_seed,
    )


def simulate_amplicon_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Run the complete amplicon read simulation pipeline.

    Args:
        config: Pipeline configuration dictionary containing tools, amplicon_params,
            pacbio_params, and read_simulation sections.
        input_fa: Path to input FASTA (haploid or diploid reference).
        human_reference: Optional path to human reference for final alignment.
            If None, returns FASTQ instead of aligned BAM.
        source_tracker: Read source tracker — not yet supported for amplicon mode.
            Passing a non-None value raises RuntimeError.
        output_config: Optional OutputConfig controlling output file placement.
            If None, output is placed alongside the input FASTA.

    Returns:
        Path to final output file (aligned BAM if human_reference given,
        otherwise FASTQ).

    Raises:
        RuntimeError: If source_tracker is provided, or if any pipeline stage fails.
    """
    # Reject read-source tracking for v1
    if source_tracker is not None:
        raise ReadSimulationError(
            "Read source tracking is not yet supported for amplicon simulation. "
            "Remove --track-read-source to proceed."
        )

    # Extract configuration
    from typing import cast

    from ..type_defs import AmpliconConfig, PacbioConfig, ReadSimulationConfig

    _amplicon_raw = config.get("amplicon_params", {})
    _pacbio_raw = config.get("pacbio_params", {})
    _rs_raw = config.get("read_simulation", {})
    amplicon_params = cast(AmpliconConfig, _amplicon_raw)
    pacbio_params = cast(PacbioConfig, _pacbio_raw)
    tools = config.get("tools", {})
    rs_config = cast(ReadSimulationConfig, _rs_raw)

    # Tool paths
    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    ccs_cmd = tools.get("ccs", "ccs")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    # PacBio params
    model_type = pacbio_params["model_type"]
    model_file = pacbio_params["model_file"]
    pass_num = pacbio_params.get("pass_num", 3)
    min_passes = pacbio_params.get("min_passes", 3)
    min_rq = pacbio_params.get("min_rq", 0.99)
    threads = pacbio_params.get("threads", 4)
    seed = pacbio_params.get("seed")
    accuracy_mean = pacbio_params.get("accuracy_mean", 0.85)

    # Amplicon params
    forward_primer = amplicon_params["forward_primer"]
    reverse_primer = amplicon_params["reverse_primer"]
    expected_range = amplicon_params.get("expected_product_range")
    expected_range_tuple: tuple[int, int] | None = None
    if expected_range is not None:
        if len(expected_range) != 2:
            raise ValueError(
                "amplicon_params.expected_product_range must contain exactly "
                "2 values: [min_product_size, max_product_size]"
            )
        expected_range_tuple = (int(expected_range[0]), int(expected_range[1]))

    # Coverage
    # Prefer read_simulation coverage, fall back to pacbio_params, then default 30.
    # Use explicit None checks to preserve coverage=0 if ever set intentionally.
    _rs_cov = rs_config.get("coverage")
    _pb_cov = pacbio_params.get("coverage")
    total_coverage: int = int(
        _rs_cov if _rs_cov is not None else (_pb_cov if _pb_cov is not None else 30)
    )

    # PCR bias config
    pcr_bias_config = amplicon_params.get("pcr_bias", {})

    # Output paths
    output_dir, output_base = resolve_pipeline_outputs(input_fa, _rs_raw, output_config)

    start_time = datetime.now()
    intermediate_files: list[str] = []

    try:
        with tempfile.TemporaryDirectory(prefix="amplicon_sim_") as temp_dir:
            temp_path = Path(temp_dir)

            logging.info("=" * 80)
            logging.info("AMPLICON SIMULATION PIPELINE")
            logging.info("=" * 80)

            # STAGES 1-4: Shared amplicon extraction and preparation
            prep = extract_and_prepare_amplicons(
                input_fa=input_fa,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                total_coverage=total_coverage,
                work_dir=temp_path,
                expected_product_range=expected_range_tuple,
                pcr_bias_config=pcr_bias_config,
                seed=seed,
            )

            intermediate_files.extend(prep.intermediate_files)
            template_fastas = [str(t) for t in prep.allele_templates]

            # STAGE 5: PBSIM3 template mode simulation
            logging.info("STAGE 5: Running PBSIM3 template mode simulation")

            clr_bam_groups: list[list[str]] = []
            for i, template_fa in enumerate(template_fastas, 1):
                prefix = str(temp_path / f"clr_hap{i}")
                hap_seed = (seed + i) if seed is not None else None

                bams = run_pbsim3_template_simulation(
                    pbsim3_cmd=pbsim3_cmd,
                    samtools_cmd=samtools_cmd,
                    template_fasta=template_fa,
                    model_type=model_type,
                    model_file=model_file,
                    output_prefix=prefix,
                    pass_num=pass_num,
                    accuracy_mean=accuracy_mean,
                    seed=hap_seed,
                )
                clr_bam_groups.append(bams)
                intermediate_files.extend(bams)
                logging.info("  Haplotype %d: %d CLR BAMs", i, len(bams))

            # STAGE 6: CCS consensus generation
            logging.info("STAGE 6: Generating HiFi consensus with CCS")

            hifi_bams = []
            for i, clr_bams in enumerate(clr_bam_groups, 1):
                for j, clr_bam in enumerate(clr_bams, 1):
                    hifi_out = str(temp_path / f"hifi_hap{i}_{j:04d}.bam")
                    hifi_bam = _simulate_haplotype_amplicon(
                        clr_bam=clr_bam,
                        ccs_cmd=ccs_cmd,
                        output_path=hifi_out,
                        min_passes=min_passes,
                        min_rq=min_rq,
                        threads=threads,
                        seed=seed,
                        hap_idx=i,
                        bam_idx=j,
                    )
                    hifi_bams.append(hifi_bam)
                    intermediate_files.append(hifi_bam)

            # STAGE 7: Merge HiFi BAMs
            hifi_merged = str(output_dir / f"{output_base}_amplicon_hifi.bam")

            if len(hifi_bams) > 1:
                logging.info("STAGE 7: Merging %d HiFi BAMs", len(hifi_bams))
                hifi_merged = merge_bam_files(
                    samtools_cmd=samtools_cmd,
                    input_bams=hifi_bams,
                    output_bam=hifi_merged,
                    threads=threads,
                )
            else:
                shutil.copy(hifi_bams[0], hifi_merged)

            intermediate_files.append(hifi_merged)

            # Convert to FASTQ
            hifi_fastq = str(output_dir / f"{output_base}_amplicon_hifi.fastq")
            hifi_fastq = convert_bam_to_fastq(
                samtools_cmd=samtools_cmd,
                input_bam=hifi_merged,
                output_fastq=hifi_fastq,
                threads=threads,
            )

            # STAGE 8: Alignment (optional)
            if human_reference is None:
                logging.info("Amplicon simulation complete (no alignment)")
                logging.info("Final output: %s", hifi_fastq)
                cleanup_intermediates(intermediate_files)
                final_output = hifi_fastq
            else:
                intermediate_files.append(hifi_fastq)

                logging.info("STAGE 8: Aligning with minimap2 (map-hifi)")
                aligned_bam = str(output_dir / f"{output_base}_amplicon_aligned.bam")
                aligned_bam = align_reads_with_minimap2(
                    minimap2_cmd=minimap2_cmd,
                    samtools_cmd=samtools_cmd,
                    reference=human_reference,
                    reads_fastq=hifi_fastq,
                    output_bam=aligned_bam,
                    preset=MINIMAP2_PRESET_PACBIO_HIFI,
                    threads=threads,
                )

                cleanup_intermediates(intermediate_files)
                final_output = aligned_bam

            # Write metadata and log completion for both paths
            end_time = datetime.now()
            duration = end_time - start_time
            logging.info("=" * 80)
            logging.info(
                "Amplicon simulation complete (duration: %s)",
                str(duration).split(".")[0],
            )
            logging.info("Final output: %s", final_output)
            logging.info("=" * 80)

            create_pipeline_metadata(
                output_dir=output_dir,
                output_base=output_base,
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="PacBio",
                tools_used=["pbsim3", "ccs", "minimap2", "samtools"],
            )

            return final_output

    except (ExternalToolError, FileOperationError) as e:
        logging.error("Amplicon pipeline failed: %s", e)
        raise RuntimeError(
            f"Amplicon simulation pipeline failed.\nError: {e}\n\n"
            f"Troubleshooting:\n"
            f"  1. Verify tools installed: pbsim, ccs, samtools, minimap2\n"
            f"  2. Check model file: {model_file}\n"
            f"  3. Ensure primers bind in the reference"
        ) from e
    except Exception as e:
        logging.error("Unexpected error in amplicon pipeline: %s", e)
        raise RuntimeError(f"Amplicon pipeline failed: {e}") from e
    finally:
        try:
            if intermediate_files:
                cleanup_intermediates(intermediate_files)
        except Exception as cleanup_err:
            logging.warning("Cleanup failed (non-fatal): %s", cleanup_err)

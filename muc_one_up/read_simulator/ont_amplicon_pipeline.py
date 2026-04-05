"""ONT amplicon read simulation pipeline.

Single-pass amplicon simulation using pbsim3 template mode with ONT
error models. Shares extraction/PCR-bias stages with the PacBio
amplicon pipeline via amplicon_common.py.

Pipeline stages:
1-3. Amplicon extraction and preparation (shared)
4.   PBSIM3 template mode, ONT model, pass_num=1
5.   BAM → FASTQ conversion (no CCS — ONT is single-pass)
6.   Merge per-allele FASTQs (diploid)
7.   Align to reference (minimap2 map-ont)
8.   Write metadata, cleanup
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
from .constants import MINIMAP2_PRESET_ONT
from .pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    resolve_pipeline_outputs,
)
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_wrapper import convert_bam_to_fastq


def simulate_ont_amplicon_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Run ONT amplicon simulation pipeline.

    Uses pbsim3 template mode with ONT error models (single-pass,
    no CCS consensus). See module docstring for pipeline stages.

    Args:
        config: Pipeline config with tools, amplicon_params, pacbio_params,
            read_simulation sections.
        input_fa: Input FASTA (haploid or diploid).
        human_reference: Optional reference for alignment.
        source_tracker: Not supported — raises if non-None.
        output_config: Optional output path control.

    Returns:
        Path to final output (aligned BAM or FASTQ).
    """
    if source_tracker is not None:
        raise ReadSimulationError(
            "Read source tracking is not yet supported for amplicon simulation. "
            "Remove --track-read-source to proceed."
        )

    from typing import cast

    from ..type_defs import AmpliconConfig, PacbioConfig, ReadSimulationConfig

    amplicon_params = cast(AmpliconConfig, config.get("amplicon_params", {}))
    pacbio_params = cast(PacbioConfig, config.get("pacbio_params", {}))
    tools = config.get("tools", {})
    rs_raw = config.get("read_simulation", {})
    rs_config = cast(ReadSimulationConfig, rs_raw)

    # Tool paths (no CCS needed for ONT)
    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    # Params from config
    model_type = pacbio_params["model_type"]
    model_file = pacbio_params["model_file"]
    threads = pacbio_params.get("threads", 4)
    seed = pacbio_params.get("seed")
    accuracy_mean = pacbio_params.get("accuracy_mean", 0.85)

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

    _rs_cov = rs_config.get("coverage")
    _pb_cov = pacbio_params.get("coverage")
    total_coverage: int = int(
        _rs_cov if _rs_cov is not None else (_pb_cov if _pb_cov is not None else 30)
    )

    pcr_bias_config = amplicon_params.get("pcr_bias", {})

    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs_raw, output_config)

    start_time = datetime.now()
    intermediate_files: list[str] = []

    try:
        with tempfile.TemporaryDirectory(prefix="ont_amplicon_sim_") as temp_dir:
            temp_path = Path(temp_dir)

            logging.info("=" * 80)
            logging.info("ONT AMPLICON SIMULATION PIPELINE")
            logging.info("=" * 80)

            # STAGES 1-3: Shared amplicon extraction and preparation
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

            # STAGE 4: PBSIM3 template mode — ONT, single-pass
            logging.info("STAGE 4: Running PBSIM3 template mode (ONT, pass_num=1)")

            allele_bams: list[list[str]] = []
            for i, template_fa in enumerate(prep.allele_templates, 1):
                prefix = str(temp_path / f"ont_hap{i}")
                hap_seed = (seed + i) if seed is not None else None

                bams = run_pbsim3_template_simulation(
                    pbsim3_cmd=pbsim3_cmd,
                    samtools_cmd=samtools_cmd,
                    template_fasta=str(template_fa),
                    model_type=model_type,
                    model_file=model_file,
                    output_prefix=prefix,
                    pass_num=1,
                    accuracy_mean=accuracy_mean,
                    seed=hap_seed,
                )
                allele_bams.append(bams)
                intermediate_files.extend(bams)
                logging.info("  Haplotype %d: %d BAMs", i, len(bams))

            # STAGE 5: BAM → FASTQ conversion (no CCS for ONT)
            logging.info("STAGE 5: Converting BAMs to FASTQ (no CCS)")

            allele_fastqs: list[str] = []
            for i, bam_group in enumerate(allele_bams, 1):
                for j, bam in enumerate(bam_group, 1):
                    fq_out = str(temp_path / f"ont_hap{i}_{j:04d}.fastq")
                    convert_bam_to_fastq(
                        samtools_cmd=samtools_cmd,
                        input_bam=bam,
                        output_fastq=fq_out,
                        threads=threads,
                    )
                    allele_fastqs.append(fq_out)

            # STAGE 6: Merge FASTQs
            merged_fastq = str(output_dir / f"{output_base}_amplicon_ont.fastq")
            logging.info("STAGE 6: Merging %d FASTQs", len(allele_fastqs))

            with open(merged_fastq, "w") as outf:
                for fq in allele_fastqs:
                    with open(fq) as inf:
                        shutil.copyfileobj(inf, outf)

            # STAGE 7: Alignment (optional)
            if human_reference is None:
                logging.info("ONT amplicon simulation complete (no alignment)")
                logging.info("Final output: %s", merged_fastq)
                final_output = merged_fastq
            else:
                logging.info("STAGE 7: Aligning with minimap2 (map-ont)")
                aligned_bam = str(output_dir / f"{output_base}_amplicon_ont.bam")
                aligned_bam = align_reads_with_minimap2(
                    minimap2_cmd=minimap2_cmd,
                    samtools_cmd=samtools_cmd,
                    reference=human_reference,
                    reads_fastq=merged_fastq,
                    output_bam=aligned_bam,
                    preset=MINIMAP2_PRESET_ONT,
                    threads=threads,
                )
                final_output = aligned_bam

            # STAGE 8: Metadata
            end_time = datetime.now()
            duration = end_time - start_time
            logging.info("=" * 80)
            logging.info(
                "ONT amplicon simulation complete (duration: %s)",
                str(duration).split(".")[0],
            )
            logging.info("Final output: %s", final_output)
            logging.info("=" * 80)

            config.setdefault("read_simulation", {})["assay_type"] = "amplicon"
            create_pipeline_metadata(
                output_dir=output_dir,
                output_base=f"{output_base}_amplicon_ont",
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="ONT",
                tools_used=["pbsim3", "minimap2", "samtools"],
            )

            return final_output

    except (ExternalToolError, FileOperationError) as e:
        logging.error("ONT amplicon pipeline failed: %s", e)
        raise RuntimeError(
            f"ONT amplicon simulation pipeline failed.\nError: {e}\n\n"
            f"Troubleshooting:\n"
            f"  1. Verify tools installed: pbsim3, samtools, minimap2\n"
            f"  2. Check model file: {model_file}\n"
            f"  3. Ensure primers bind in the reference"
        ) from e
    except Exception as e:
        logging.error("Unexpected error in ONT amplicon pipeline: %s", e)
        raise RuntimeError(f"ONT amplicon pipeline failed: {e}") from e
    finally:
        try:
            if intermediate_files:
                cleanup_intermediates(intermediate_files)
        except Exception as cleanup_err:
            logging.warning("Cleanup failed (non-fatal): %s", cleanup_err)

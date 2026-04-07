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
        config: Pipeline config with tools, amplicon_params, ont_amplicon_params,
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

    from ..type_defs import AmpliconConfig, OntAmpliconConfig, ReadSimulationConfig

    amplicon_params = cast(AmpliconConfig, config.get("amplicon_params", {}))
    ont_params = cast(OntAmpliconConfig, config.get("ont_amplicon_params", {}))
    tools = config.get("tools", {})
    rs_raw = config.get("read_simulation", {})
    rs_config = cast(ReadSimulationConfig, rs_raw)

    # Tool paths (no CCS needed for ONT)
    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    # ONT-specific params — uses ont_amplicon_params (NOT pacbio_params)
    # Treat None as missing (possible from partial configs with explicit nulls)
    model_type = ont_params.get("model_type") or "errhmm"
    model_file = ont_params.get("model_file") or "reference/pbsim3/ERRHMM-ONT.model"
    threads = ont_params.get("threads") or 8
    seed = ont_params.get("seed")  # None is valid here (random seed)
    accuracy_mean = ont_params.get("accuracy_mean") or 0.95

    # Validate model file looks like an ONT model
    model_name = Path(model_file).name.upper()
    if "ONT" not in model_name and "NANOPORE" not in model_name:
        logging.warning(
            "ONT amplicon mode is using model file '%s' which does not appear "
            "to be an ONT model. The default is ERRHMM-ONT.model. "
            "Check ont_amplicon_params.model_file in config or use --model-file.",
            model_file,
        )

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
    total_coverage: int = int(_rs_cov if _rs_cov is not None else 30)

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

            allele_outputs: list[list[str]] = []
            for i, template_fa in enumerate(prep.allele_templates, 1):
                prefix = str(temp_path / f"ont_hap{i}")
                hap_seed = (seed + i) if seed is not None else None

                outputs = run_pbsim3_template_simulation(
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
                allele_outputs.append(outputs)
                intermediate_files.extend(outputs)
                logging.info("  Haplotype %d: %d output file(s)", i, len(outputs))

            # STAGE 5: Collect FASTQ reads
            # pbsim3 with pass_num=1 outputs FASTQ directly (.fq.gz);
            # with pass_num>=2 it outputs BAM (needs conversion).
            logging.info("STAGE 5: Collecting reads as FASTQ")

            allele_fastqs: list[str] = []
            for i, output_group in enumerate(allele_outputs, 1):
                for j, output_file in enumerate(output_group, 1):
                    if output_file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
                        allele_fastqs.append(output_file)
                    else:
                        fq_out = str(temp_path / f"ont_hap{i}_{j:04d}.fastq")
                        convert_bam_to_fastq(
                            samtools_cmd=samtools_cmd,
                            input_bam=output_file,
                            output_fastq=fq_out,
                            threads=threads,
                        )
                        allele_fastqs.append(fq_out)

            # STAGE 6: Merge FASTQs
            merged_fastq = str(output_dir / f"{output_base}_amplicon_ont.fastq")
            logging.info("STAGE 6: Merging %d FASTQs", len(allele_fastqs))

            from .utils.fastq_utils import merge_fastq_files

            merge_fastq_files(allele_fastqs, merged_fastq, validate_inputs=False)

            # STAGE 7: Alignment (optional)
            if human_reference is None:
                logging.info("ONT amplicon simulation complete (no alignment)")
                logging.info("Final output: %s", merged_fastq)
                final_output = merged_fastq
            else:
                intermediate_files.append(merged_fastq)
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
        keep = config.get("read_simulation", {}).get("keep_intermediate_files", False)
        if keep:
            logging.info("Keeping intermediate files (keep_intermediate_files=true)")
        else:
            try:
                if intermediate_files:
                    cleanup_intermediates(intermediate_files)
            except Exception as cleanup_err:
                logging.warning("Cleanup failed (non-fatal): %s", cleanup_err)

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
from .constants import MINIMAP2_PRESET_PACBIO_HIFI
from .pcr_bias import PCRBiasModel
from .utils import cleanup_files, write_metadata_file
from .utils.amplicon_extractor import AmpliconExtractor
from .utils.reference_utils import extract_haplotypes, is_diploid_reference
from .utils.template_generator import generate_template_fasta
from .wrappers.ccs_wrapper import run_ccs_consensus
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_wrapper import convert_bam_to_fastq, merge_bam_files


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
    amplicon_params = config.get("amplicon_params", {})
    pacbio_params = config.get("pacbio_params", {})
    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})

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
    expected_range_tuple = tuple(expected_range) if expected_range else None

    # Coverage
    total_coverage = rs_config.get("coverage", pacbio_params.get("coverage", 30))

    # PCR bias model
    pcr_bias_config = amplicon_params.get("pcr_bias", {})
    pcr_model = PCRBiasModel.from_config(pcr_bias_config)

    # Output paths
    input_path = Path(input_fa)
    if output_config is not None:
        output_dir = output_config.out_dir
        output_base = output_config.out_base
    else:
        output_dir = input_path.parent
        output_base = input_path.stem

    output_dir.mkdir(parents=True, exist_ok=True)

    start_time = datetime.now()
    intermediate_files: list[str] = []

    try:
        with tempfile.TemporaryDirectory(prefix="amplicon_sim_") as temp_dir:
            temp_path = Path(temp_dir)

            logging.info("=" * 80)
            logging.info("AMPLICON SIMULATION PIPELINE")
            logging.info("=" * 80)

            # STAGE 1: Haplotype detection and extraction
            diploid = is_diploid_reference(input_fa)

            if diploid:
                logging.info("STAGE 1: Extracting haplotypes from diploid reference")
                hap1_fa, hap2_fa = extract_haplotypes(input_fa, temp_path, base_name="amplicon_sim")
                haplotype_fastas = [hap1_fa, hap2_fa]
            else:
                logging.info("STAGE 1: Haploid reference — skipping extraction")
                haplotype_fastas = [input_fa]

            # STAGE 2: Amplicon extraction per haplotype
            logging.info("STAGE 2: Extracting amplicons via primer binding sites")

            extractor = AmpliconExtractor(
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                expected_product_range=expected_range_tuple,
            )

            amplicon_results = []
            for i, hap_fa in enumerate(haplotype_fastas, 1):
                amp_out = str(temp_path / f"amplicon_hap{i}.fa")
                result = extractor.extract(hap_fa, amp_out)
                amplicon_results.append(result)
                logging.info("  Haplotype %d amplicon: %d bp", i, result.length)

            # STAGE 3: PCR bias coverage split
            if diploid:
                logging.info("STAGE 3: Computing PCR bias coverage split")
                n1, n2 = pcr_model.compute_coverage_split(
                    total_coverage,
                    amplicon_results[0].length,
                    amplicon_results[1].length,
                    seed=seed,
                )
                allele_counts = [n1, n2]
                logging.info(
                    "  Coverage split: allele1=%d, allele2=%d (total=%d)",
                    n1,
                    n2,
                    total_coverage,
                )
            else:
                logging.info("STAGE 3: Haploid — full coverage to single allele")
                allele_counts = [total_coverage]

            # STAGE 4: Template FASTA generation
            logging.info("STAGE 4: Generating template FASTAs")

            template_fastas = []
            for i, (amp_result, count) in enumerate(
                zip(amplicon_results, allele_counts, strict=False), 1
            ):
                template_out = str(temp_path / f"template_hap{i}.fa")
                generate_template_fasta(amp_result.fasta_path, count, template_out)
                template_fastas.append(template_out)
                logging.info(
                    "  Haplotype %d: %d copies of %d bp amplicon",
                    i,
                    count,
                    amp_result.length,
                )

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
                    hap_seed = (seed + i * 100 + j) if seed is not None else None
                    hifi_bam = run_ccs_consensus(
                        ccs_cmd=ccs_cmd,
                        input_bam=clr_bam,
                        output_bam=hifi_out,
                        min_passes=min_passes,
                        min_rq=min_rq,
                        threads=threads,
                        seed=hap_seed,
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
                cleanup_files(intermediate_files)
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

                cleanup_files(intermediate_files)
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

            write_metadata_file(
                output_dir=str(output_dir),
                output_base=output_base,
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="PacBio-Amplicon",
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
                cleanup_files(intermediate_files)
        except Exception as cleanup_err:
            logging.warning("Cleanup failed (non-fatal): %s", cleanup_err)

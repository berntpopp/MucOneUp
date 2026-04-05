"""Shared amplicon extraction and preparation stages.

Used by both PacBio and ONT amplicon pipelines. Handles:
- Diploid detection and haplotype extraction
- Primer-based amplicon extraction
- PCR bias coverage split
- Template FASTA generation
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .pcr_bias import PCRBiasModel
from .utils.amplicon_extractor import AmpliconExtractor
from .utils.reference_utils import extract_haplotypes, is_diploid_reference
from .utils.template_generator import generate_template_fasta


@dataclass
class AmpliconPrep:
    """Result of amplicon extraction and preparation stages."""

    allele_templates: list[Path]
    allele_coverages: list[int]
    output_dir: Path
    output_base: str
    intermediate_files: list[str] = field(default_factory=list)
    is_diploid: bool = False


def extract_and_prepare_amplicons(
    *,
    input_fa: str,
    forward_primer: str,
    reverse_primer: str,
    total_coverage: int,
    work_dir: Path,
    expected_product_range: tuple[int, int] | None = None,
    pcr_bias_config: dict[str, Any] | None = None,
    seed: int | None = None,
) -> AmpliconPrep:
    """Extract amplicons and prepare template FASTAs for simulation.

    Stages:
    1. Detect diploid, split haplotypes if needed
    2. Extract amplicon via primer binding sites
    3. Compute PCR bias coverage split (diploid only)
    4. Generate template FASTAs (N copies per allele)

    Args:
        input_fa: Path to input FASTA (haploid or diploid).
        forward_primer: Forward primer sequence.
        reverse_primer: Reverse primer sequence.
        total_coverage: Total number of template molecules.
        work_dir: Working directory for intermediate files.
        expected_product_range: Optional (min, max) expected amplicon size.
        pcr_bias_config: PCR bias config dict (preset, stochastic, etc.).
        seed: Random seed for reproducibility.

    Returns:
        AmpliconPrep with template paths, coverage per allele, etc.
    """
    intermediate_files: list[str] = []

    # STAGE 1: Haplotype detection and extraction
    diploid = is_diploid_reference(input_fa)

    if diploid:
        logging.info("STAGE 1: Extracting haplotypes from diploid reference")
        hap1_fa, hap2_fa = extract_haplotypes(input_fa, work_dir, base_name="amplicon_sim")
        haplotype_fastas = [hap1_fa, hap2_fa]
    else:
        logging.info("STAGE 1: Haploid reference — skipping extraction")
        haplotype_fastas = [input_fa]

    # STAGE 2: Amplicon extraction per haplotype
    logging.info("STAGE 2: Extracting amplicons via primer binding sites")

    extractor = AmpliconExtractor(
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        expected_product_range=expected_product_range,
    )

    amplicon_results = []
    for i, hap_fa in enumerate(haplotype_fastas, 1):
        amp_out = str(work_dir / f"amplicon_hap{i}.fa")
        result = extractor.extract(hap_fa, amp_out)
        amplicon_results.append(result)
        logging.info("  Haplotype %d amplicon: %d bp", i, result.length)

    # STAGE 3: PCR bias coverage split
    pcr_model = PCRBiasModel.from_config(pcr_bias_config or {})

    if diploid:
        logging.info("STAGE 3: Computing PCR bias coverage split")
        n1, n2 = pcr_model.compute_coverage_split(
            total_coverage,
            amplicon_results[0].length,
            amplicon_results[1].length,
            seed=seed,
        )
        allele_counts = [max(1, n1), max(1, n2)]
        if n1 == 0 or n2 == 0:
            logging.warning(
                "PCR bias produced 0 reads for one allele at coverage=%d. "
                "Each allele will receive at least 1 read.",
                total_coverage,
            )
        logging.info(
            "  Coverage split: allele1=%d, allele2=%d (total=%d)",
            allele_counts[0],
            allele_counts[1],
            total_coverage,
        )
    else:
        logging.info("STAGE 3: Haploid — full coverage to single allele")
        allele_counts = [total_coverage]

    # STAGE 4: Template FASTA generation
    logging.info("STAGE 4: Generating template FASTAs")

    template_fastas: list[Path] = []
    for i, (amp_result, count) in enumerate(zip(amplicon_results, allele_counts, strict=False), 1):
        template_out = work_dir / f"template_hap{i}.fa"
        generate_template_fasta(amp_result.fasta_path, count, str(template_out))
        template_fastas.append(template_out)
        logging.info(
            "  Haplotype %d: %d copies of %d bp amplicon",
            i,
            count,
            amp_result.length,
        )

    return AmpliconPrep(
        allele_templates=template_fastas,
        allele_coverages=allele_counts,
        output_dir=work_dir,
        output_base="amplicon",
        intermediate_files=intermediate_files,
        is_diploid=diploid,
    )

"""Genomic/targeted ONT simulation from sampled fragments with per-read truth.

Alternative to the NanoSim ``ont`` simulator. Fragments are drawn from the
simulated haplotypes (optionally extended with extra flank sequence), passed
through the same molecule model (strand, homopolymer stutter) and pbsim3
template path as amplicon simulation, so every read carries exact truth.

Stages:
1. Load haplotypes (+ optional ``flank_fasta`` records ``left``/``right``)
2. Sample fragments: uniform haplotype and start, log-normal length
3. pbsim3 template mode (ONT error model), relabel reads, write truth
4. Optional alignment (minimap2 map-ont), metadata
"""

from __future__ import annotations

import logging
import random
import tempfile
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

from Bio import SeqIO

from .constants import MINIMAP2_PRESET_ONT
from .molecule_pipeline import PbsimRun, sequencer_for, simulate_molecule_reads
from .molecules import MoleculeModel, build_fragment_molecules
from .pipeline_utils import create_pipeline_metadata, resolve_pipeline_outputs
from .read_profiles import FragmentModel, active_read_profile
from .wrappers.minimap2_wrapper import align_reads_with_minimap2

if TYPE_CHECKING:
    from .output_config import OutputConfig

DEFAULT_ONT_MODEL = "reference/pbsim3/QSHMM-ONT-HQ.model"


def load_fragment_sources(input_fa: str, flank_fasta: str | None) -> list[str]:
    """Haplotype sequences, each wrapped in the optional ``left``/``right`` flank records."""
    haplotypes = [str(rec.seq).upper() for rec in SeqIO.parse(input_fa, "fasta")]
    if not haplotypes:
        raise ValueError(f"No sequences in {input_fa}")
    left = right = ""
    if flank_fasta:
        flanks = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(flank_fasta, "fasta")}
        unknown = set(flanks) - {"left", "right"}
        if unknown:
            raise ValueError(
                f"flank_fasta records must be named 'left'/'right', got {sorted(unknown)}"
            )
        left, right = flanks.get("left", ""), flanks.get("right", "")
    return [left + hap + right for hap in haplotypes]


def _fragment_settings(
    config: dict[str, Any],
) -> tuple[MoleculeModel, FragmentModel, dict[str, Any]]:
    profile = active_read_profile(config)
    molecules = profile.molecules if profile else MoleculeModel(forward_frac=0.5)
    fragments = profile.fragments if profile else FragmentModel()
    params = dict(config.get("ont_fragment_params", {}))
    if "length_median" in params or "length_sigma" in params:
        fragments = FragmentModel(
            params.get("length_median", fragments.length_median),
            params.get("length_sigma", fragments.length_sigma),
        )
    return molecules, fragments, params


def _n_reads(params: dict[str, Any], rs: dict[str, Any], sources: list[str], median: float) -> int:
    if params.get("n_reads"):
        return int(params["n_reads"])
    coverage = float(rs.get("coverage") or 30)
    mean_source = sum(len(s) for s in sources) / len(sources)
    return max(1, round(coverage * mean_source * len(sources) / median))


def simulate_ont_fragment_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Simulate truth-tracked genomic ONT reads; returns the FASTQ or aligned BAM path."""
    del source_tracker  # truth is always written by this simulator
    rs = config.get("read_simulation", {})
    tools = config.get("tools", {})
    ont = config.get("ont_amplicon_params", {})  # pbsim3 ONT settings shared with amplicon mode
    molecules_model, fragments, params = _fragment_settings(config)
    seed = params.get("seed", ont.get("seed"))
    sources = load_fragment_sources(input_fa, params.get("flank_fasta"))
    n_reads = _n_reads(params, rs, sources, fragments.length_median)
    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs, output_config)
    fastq = output_dir / f"{output_base}_ont_fragments.fastq"
    truth = output_dir / f"{output_base}_read_truth.tsv.gz"
    run = PbsimRun(
        pbsim3_cmd=tools.get("pbsim3", "pbsim"),
        samtools_cmd=tools.get("samtools", "samtools"),
        model_type=ont.get("model_type") or "qshmm",
        model_file=ont.get("model_file") or DEFAULT_ONT_MODEL,
        accuracy_mean=ont.get("accuracy_mean") or 0.95,
        difference_ratio=ont.get("difference_ratio"),
        threads=ont.get("threads") or 4,
    )
    start = datetime.now()
    logging.info(
        "ONT fragment simulation: %d reads, median %.0f bp", n_reads, fragments.length_median
    )
    molecules = build_fragment_molecules(
        sources,
        n_reads,
        fragments.length_median,
        fragments.length_sigma,
        molecules_model,
        random.Random(seed),
    )
    with tempfile.TemporaryDirectory(prefix="ont_fragment_sim_") as tmp:
        simulate_molecule_reads(
            molecules, sequencer_for(config, run), Path(tmp), fastq, truth, output_base, seed
        )
    final_output = str(fastq)
    if human_reference:
        final_output = align_reads_with_minimap2(
            minimap2_cmd=tools.get("minimap2", "minimap2"),
            samtools_cmd=run.samtools_cmd,
            reference=human_reference,
            reads_fastq=str(fastq),
            output_bam=str(output_dir / f"{output_base}_ont_fragments.bam"),
            preset=MINIMAP2_PRESET_ONT,
            threads=run.threads,
        )
    create_pipeline_metadata(
        output_dir=output_dir,
        output_base=f"{output_base}_ont_fragments",
        config=config,
        start_time=start,
        end_time=datetime.now(),
        platform="ONT",
        tools_used=["pbsim3", "minimap2", "samtools"],
    )
    return final_output

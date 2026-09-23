"""Simulate truth-tracked reads from template molecules with pbsim3 (and ccs for HiFi).

One pbsim3 template-mode run covers all molecules. Single-pass output (ONT) is
FASTQ; multi-pass output (PacBio) is turned into HiFi reads with ccs. Reads are
then relabelled from the pbsim3 MAF so every read carries molecule truth.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from .molecules import Molecule
from .read_truth import parse_maf_read_templates, relabel_reads, template_id
from .wrappers.ccs_wrapper import run_ccs_consensus
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_convert import convert_bam_to_fastq

_READ_PREFIX = "mol"


@dataclass(frozen=True)
class PbsimRun:
    """Tool commands and error-model settings for one molecule simulation."""

    pbsim3_cmd: str
    samtools_cmd: str
    model_type: str
    model_file: str
    pass_num: int = 1
    accuracy_mean: float = 0.95
    accuracy_sd: float | None = None
    difference_ratio: str | None = None
    ccs_cmd: str | None = None
    min_passes: int = 3
    min_rq: float = 0.99
    threads: int = 4

    def __post_init__(self) -> None:
        if self.pass_num > 1 and not self.ccs_cmd:
            raise ValueError("multi-pass (HiFi) simulation requires ccs_cmd")


def write_molecule_templates(molecules: Sequence[Molecule], path: Path) -> Path:
    """Write one FASTA record per molecule, named by molecule id."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for molecule in molecules:
            handle.write(f">{template_id(molecule.id)}\n{molecule.seq}\n")
    return path


def _to_fastqs(outputs: list[str], run: PbsimRun, work_dir: Path, seed: int | None) -> list[Path]:
    if run.pass_num == 1:
        return [Path(p) for p in outputs]
    assert run.ccs_cmd is not None  # guaranteed by PbsimRun.__post_init__
    fastqs = []
    for index, clr_bam in enumerate(outputs, start=1):
        hifi_bam = run_ccs_consensus(
            ccs_cmd=run.ccs_cmd,
            input_bam=clr_bam,
            output_bam=str(work_dir / f"hifi_{index:04d}.bam"),
            min_passes=run.min_passes,
            min_rq=run.min_rq,
            threads=run.threads,
            seed=None if seed is None else seed + index,
        )
        fastqs.append(
            Path(
                convert_bam_to_fastq(
                    samtools_cmd=run.samtools_cmd,
                    input_bam=hifi_bam,
                    output_fastq=str(work_dir / f"hifi_{index:04d}.fastq"),
                    threads=run.threads,
                )
            )
        )
    return fastqs


def simulate_molecule_reads(
    molecules: Sequence[Molecule],
    run: PbsimRun,
    work_dir: Path,
    out_fastq: Path,
    truth_tsv: Path,
    base: str,
    seed: int | None,
) -> int:
    """Simulate reads for ``molecules``; write relabelled FASTQ and truth. Returns read count."""
    work_dir.mkdir(parents=True, exist_ok=True)
    templates = write_molecule_templates(molecules, work_dir / "molecules.fa")
    prefix = work_dir / "molecules"
    outputs = run_pbsim3_template_simulation(
        pbsim3_cmd=run.pbsim3_cmd,
        samtools_cmd=run.samtools_cmd,
        template_fasta=str(templates),
        model_type=run.model_type,
        model_file=run.model_file,
        output_prefix=str(prefix),
        pass_num=run.pass_num,
        accuracy_mean=run.accuracy_mean,
        seed=seed,
        accuracy_sd=run.accuracy_sd,
        difference_ratio=run.difference_ratio,
        id_prefix=_READ_PREFIX,
    )
    mapping: dict[str, str] = {}
    for maf in sorted(work_dir.glob(f"{prefix.name}*.maf.gz")):
        mapping.update(parse_maf_read_templates(maf))
    fastqs = _to_fastqs(outputs, run, work_dir, seed)
    count = relabel_reads(fastqs, mapping, {m.id: m for m in molecules}, base, out_fastq, truth_tsv)
    logging.info("Simulated %d truth-tracked reads from %d molecules", count, len(molecules))
    return count

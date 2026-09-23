"""Simulate truth-tracked reads from template molecules.

A *sequencer* turns molecules into reads and reports which template each read
came from (Strategy pattern):

* :class:`PbsimRun` - one pbsim3 template-mode run for all molecules; single-pass
  output (ONT) is FASTQ, multi-pass output (PacBio) becomes HiFi via ccs; the
  read -> template map comes from the pbsim3 MAF.
* :class:`EmpiricalSequencer` - the calibrated in-process error channel
  (:mod:`.empirical_errors`), used when a read profile defines ``errors``.

Reads are then relabelled so every read carries molecule truth.
"""

from __future__ import annotations

import hashlib
import logging
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

from ..exceptions import ReadSimulationError
from .empirical_errors import EmpiricalErrorModel, apply_errors
from .molecules import Molecule
from .read_profiles import active_read_profile
from .read_truth import parse_maf_read_templates, relabel_reads, template_id
from .wrappers.ccs_wrapper import run_ccs_consensus
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_convert import convert_bam_to_fastq

_READ_PREFIX = "mol"


class Sequencer(Protocol):
    """Turns molecules into reads; returns FASTQ paths and a read -> template map."""

    def sequence(
        self, molecules: Sequence[Molecule], work_dir: Path, seed: int | None
    ) -> tuple[list[Path], dict[str, str]]: ...


@dataclass(frozen=True)
class PbsimRun:
    """Tool commands and error-model settings for one molecule simulation."""

    pbsim3_cmd: str
    samtools_cmd: str
    model_type: str
    model_file: str
    pass_num: int = 1
    accuracy_mean: float = 0.95
    difference_ratio: str | None = None
    ccs_cmd: str | None = None
    min_passes: int = 3
    min_rq: float = 0.99
    threads: int = 4

    def __post_init__(self) -> None:
        if self.pass_num > 1 and not self.ccs_cmd:
            raise ValueError("multi-pass (HiFi) simulation requires ccs_cmd")

    def sequence(
        self, molecules: Sequence[Molecule], work_dir: Path, seed: int | None
    ) -> tuple[list[Path], dict[str, str]]:
        templates = write_molecule_templates(molecules, work_dir / "molecules.fa")
        prefix = work_dir / "molecules"
        outputs = run_pbsim3_template_simulation(
            pbsim3_cmd=self.pbsim3_cmd,
            samtools_cmd=self.samtools_cmd,
            template_fasta=str(templates),
            model_type=self.model_type,
            model_file=self.model_file,
            output_prefix=str(prefix),
            pass_num=self.pass_num,
            accuracy_mean=self.accuracy_mean,
            seed=seed,
            difference_ratio=self.difference_ratio,
            id_prefix=_READ_PREFIX,
        )
        mapping: dict[str, str] = {}
        for maf in sorted(work_dir.glob(f"{prefix.name}*.maf.gz")):
            mapping.update(parse_maf_read_templates(maf))
        return _to_fastqs(outputs, self, work_dir, seed), mapping


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
    sequencer: Sequencer,
    work_dir: Path,
    out_fastq: Path,
    truth_tsv: Path,
    base: str,
    seed: int | None,
) -> int:
    """Sequence ``molecules``; write relabelled FASTQ and truth. Returns read count."""
    work_dir.mkdir(parents=True, exist_ok=True)
    fastqs, mapping = sequencer.sequence(molecules, work_dir, seed)
    count = relabel_reads(fastqs, mapping, {m.id: m for m in molecules}, base, out_fastq, truth_tsv)
    if count == 0:
        raise ReadSimulationError(
            f"Simulation produced 0 reads from {len(molecules)} molecules; check coverage, "
            "template lengths and ccs filters (min_passes, min_rq)"
        )
    logging.info("Simulated %d truth-tracked reads from %d molecules", count, len(molecules))
    return count


def stage_seed(seed: int | None, stage: str) -> int | None:
    """Stable per-stage seed derived from the user seed (None stays unseeded).

    The molecule model uses ``random.Random(seed)``; stages that draw from their
    own Python RNG use a derived seed so their stream is not a replay of the
    molecule stream.
    """
    if seed is None:
        return None
    digest = hashlib.sha256(f"{seed}:{stage}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


@dataclass(frozen=True)
class EmpiricalSequencer:
    """Calibrated in-process error channel; one read per molecule."""

    model: EmpiricalErrorModel

    def sequence(
        self, molecules: Sequence[Molecule], work_dir: Path, seed: int | None
    ) -> tuple[list[Path], dict[str, str]]:
        rng = random.Random(stage_seed(seed, "errors"))
        fastq = work_dir / "empirical.fastq"
        with open(fastq, "w") as handle:
            for molecule in molecules:
                read, qual = apply_errors(molecule.seq, self.model, rng)
                handle.write(f"@{template_id(molecule.id)}\n{read}\n+\n{qual}\n")
        names = {template_id(m.id): template_id(m.id) for m in molecules}
        return [fastq], names


def sequencer_for(config: Mapping[str, Any], pbsim: PbsimRun) -> Sequencer:
    """Empirical channel when the active read profile defines ``errors``, else pbsim3."""
    profile = active_read_profile(config)
    if profile is not None and profile.errors is not None:
        return EmpiricalSequencer(profile.errors)
    return pbsim

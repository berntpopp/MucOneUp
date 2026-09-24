"""Tests for simulating reads from molecules (pbsim3/ccs mocked at the tool boundary)."""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.read_simulator.molecule_pipeline import (
    PbsimRun,
    simulate_molecule_reads,
    write_molecule_templates,
)
from muc_one_up.read_simulator.molecules import Molecule

MOD = "muc_one_up.read_simulator.molecule_pipeline"
MOLS = [
    Molecule(1, 1, "full", "+", "ACGTACGT", 0, 8),
    Molecule(2, 2, "full", "-", "TTGGCCAA", 0, 8),
]


def _gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_write_templates_uses_molecule_ids(tmp_path: Path) -> None:
    path = write_molecule_templates(MOLS, tmp_path / "t.fa")
    assert path.read_text() == ">m0000001\nACGTACGT\n>m0000002\nTTGGCCAA\n"


def test_ont_single_pass_reads_are_relabelled_with_truth(tmp_path: Path) -> None:
    def fake_pbsim(**kw):
        prefix = kw["output_prefix"]
        _gz(Path(f"{prefix}.fq.gz"), "@mol_2\nTTGG\n+\nIIII\n@mol_1\nACGT\n+\nIIII\n")
        _gz(
            Path(f"{prefix}.maf.gz"),
            "a\ns m0000001 0 8 + 8 A\ns mol_1 0 4 + 4 A\n\n"
            "a\ns m0000002 0 8 + 8 A\ns mol_2 0 4 + 4 A\n\n",
        )
        return [f"{prefix}.fq.gz"]

    run = PbsimRun("pbsim", "samtools", "qshmm", "m.model", difference_ratio="33:31:36")
    with patch(f"{MOD}.run_pbsim3_template_simulation", side_effect=fake_pbsim) as pbsim:
        n = simulate_molecule_reads(
            MOLS, run, tmp_path, tmp_path / "o.fastq", tmp_path / "t.tsv.gz", "s", 7
        )
    assert n == 2
    assert pbsim.call_args.kwargs["difference_ratio"] == "33:31:36"
    assert pbsim.call_args.kwargs["pass_num"] == 1
    assert (tmp_path / "o.fastq").read_text().splitlines()[0] == "@s_h2_m0000002"


def test_hifi_runs_ccs_and_converts(tmp_path: Path) -> None:
    def fake_pbsim(**kw):
        prefix = kw["output_prefix"]
        Path(f"{prefix}.bam").write_bytes(b"BAM")
        _gz(Path(f"{prefix}.maf.gz"), "a\ns m0000001 0 8 + 8 A\ns mol/1/0 0 4 + 4 A\n\n")
        return [f"{prefix}.bam"]

    def fake_convert(**kw):
        Path(kw["output_fastq"]).write_text("@mol/1/ccs\nACGT\n+\nIIII\n")
        return kw["output_fastq"]

    run = PbsimRun("pbsim", "samtools", "errhmm", "m.model", pass_num=10, ccs_cmd="ccs")
    with (
        patch(f"{MOD}.run_pbsim3_template_simulation", side_effect=fake_pbsim),
        patch(f"{MOD}.run_ccs_consensus", side_effect=lambda **kw: kw["output_bam"]) as ccs,
        patch(f"{MOD}.convert_bam_to_fastq", side_effect=fake_convert),
    ):
        n = simulate_molecule_reads(
            MOLS[:1], run, tmp_path, tmp_path / "o.fastq", tmp_path / "t.tsv.gz", "s", 1
        )
    assert n == 1 and ccs.call_count == 1


def test_multipass_requires_ccs_command() -> None:
    with pytest.raises(ValueError, match="ccs_cmd"):
        PbsimRun("pbsim", "samtools", "errhmm", "m.model", pass_num=10)


ERRORS = {
    "mismatch_rate": 0.01,
    "insertion_rate": 0.005,
    "deletion_rate": 0.005,
    "insertion_len_pmf": {"1": 1.0},
    "deletion_len_pmf": {"1": 1.0},
    "read_error_sigma": 0.0,
}


def test_empirical_sequencer_end_to_end_without_external_tools(tmp_path: Path) -> None:
    from muc_one_up.read_simulator.empirical_errors import EmpiricalErrorModel
    from muc_one_up.read_simulator.molecule_pipeline import EmpiricalSequencer

    mols = [
        Molecule(i, 1 + i % 2, "full", "+-"[i % 2], "ACGTTGCA" * 50, 0, 400) for i in range(1, 21)
    ]
    seq = EmpiricalSequencer(EmpiricalErrorModel.from_dict(ERRORS))
    n = simulate_molecule_reads(
        mols, seq, tmp_path, tmp_path / "o.fq", tmp_path / "t.tsv.gz", "e", 3
    )
    assert n == 20
    names = (tmp_path / "o.fq").read_text().splitlines()[::4]
    assert names[:2] == ["@e_h2_m0000001", "@e_h1_m0000002"]


def test_sequencer_for_selects_engine(tmp_path: Path) -> None:
    import json

    from muc_one_up.read_simulator.molecule_pipeline import EmpiricalSequencer, sequencer_for

    pbsim = PbsimRun("pbsim", "samtools", "qshmm", "m.model")
    assert sequencer_for({}, pbsim) is pbsim
    profile = tmp_path / "p.json"
    profile.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "name": "p",
                "platform": "ont",
                "errors": ERRORS,
                "molecules": {  # an errors model needs stutter plus a fallback (#132)
                    "stutter": {"C7|+": {"-1": 0.2, "0": 0.8}},
                    "stutter_fallback": {
                        "rule": "log_odds_linear",
                        "log_odds_slope_per_base": 0.5,
                        "generic": {"ref_len": 3, "pmf": {"-1": 0.05, "0": 0.95}},
                    },
                },
            }
        )
    )
    chosen = sequencer_for({"read_model": {"profile": str(profile)}}, pbsim)
    assert isinstance(chosen, EmpiricalSequencer) and chosen.model.mismatch_rate == 0.01


def test_zero_reads_is_an_error(tmp_path: Path) -> None:
    """An empty simulation must fail loudly, not exit 0 with an empty FASTQ (#113)."""
    from muc_one_up.exceptions import ReadSimulationError
    from muc_one_up.read_simulator.empirical_errors import EmpiricalErrorModel
    from muc_one_up.read_simulator.molecule_pipeline import EmpiricalSequencer

    seq = EmpiricalSequencer(EmpiricalErrorModel.from_dict(ERRORS))
    with pytest.raises(ReadSimulationError, match="0 reads from 0 molecules"):
        simulate_molecule_reads([], seq, tmp_path, tmp_path / "o.fq", tmp_path / "t.tsv.gz", "e", 3)


def test_all_reads_dropped_is_an_error(tmp_path: Path) -> None:
    """e.g. ccs rejecting every ZMW: pbsim3/ccs succeed but nothing is emitted (#113)."""
    from muc_one_up.exceptions import ReadSimulationError

    def fake_pbsim(**kw):
        prefix = kw["output_prefix"]
        _gz(Path(f"{prefix}.fq.gz"), "")
        _gz(Path(f"{prefix}.maf.gz"), "")
        return [f"{prefix}.fq.gz"]

    run = PbsimRun("pbsim", "samtools", "qshmm", "m.model")
    with (
        patch(f"{MOD}.run_pbsim3_template_simulation", side_effect=fake_pbsim),
        pytest.raises(ReadSimulationError, match="0 reads from 2 molecules"),
    ):
        simulate_molecule_reads(
            MOLS, run, tmp_path, tmp_path / "o.fq", tmp_path / "t.tsv.gz", "b", 1
        )


def test_empirical_channel_uses_a_stream_independent_of_the_molecule_rng(tmp_path: Path) -> None:
    """Molecules are built with Random(seed); errors must not replay that stream (#121)."""
    import random

    from muc_one_up.read_simulator.empirical_errors import EmpiricalErrorModel
    from muc_one_up.read_simulator.molecule_pipeline import EmpiricalSequencer, stage_seed

    seen = []

    def spy(seq, model, rng):
        seen.append(rng.random())
        return seq, "I" * len(seq)

    with patch(f"{MOD}.apply_errors", side_effect=spy):
        EmpiricalSequencer(EmpiricalErrorModel.from_dict(ERRORS)).sequence(MOLS[:1], tmp_path, 5)
    assert seen[0] != random.Random(5).random()
    assert seen[0] == random.Random(stage_seed(5, "errors")).random()
    assert stage_seed(None, "errors") is None
    assert stage_seed(5, "errors") == stage_seed(5, "errors") != stage_seed(6, "errors")

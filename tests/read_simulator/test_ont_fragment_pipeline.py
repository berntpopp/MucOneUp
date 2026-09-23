"""Tests for the pbsim3 fragment-based genomic ONT simulator (#107)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.read_simulator.ont_fragment_pipeline import (
    load_fragment_sources,
    simulate_ont_fragment_pipeline,
)
from muc_one_up.read_simulator.output_config import OutputConfig

MOD = "muc_one_up.read_simulator.ont_fragment_pipeline"


@pytest.fixture
def diploid(tmp_path: Path) -> Path:
    fa = tmp_path / "s.simulated.fa"
    fa.write_text(">haplotype_1\n" + "A" * 3000 + "\n>haplotype_2\n" + "C" * 3000 + "\n")
    return fa


def test_sources_wrap_flanks(tmp_path: Path, diploid: Path) -> None:
    flanks = tmp_path / "f.fa"
    flanks.write_text(">left\nGG\n>right\nTT\n")
    sources = load_fragment_sources(str(diploid), str(flanks))
    assert [s[:3] for s in sources] == ["GGA", "GGC"] and all(s.endswith("TT") for s in sources)


def test_flank_records_must_be_named(tmp_path: Path, diploid: Path) -> None:
    bad = tmp_path / "f.fa"
    bad.write_text(">chr1\nGG\n")
    with pytest.raises(ValueError, match="left"):
        load_fragment_sources(str(diploid), str(bad))


def _run(tmp_path: Path, diploid: Path, config: dict):
    with (
        patch(f"{MOD}.simulate_molecule_reads", return_value=0) as sim,
        patch(f"{MOD}.create_pipeline_metadata"),
        patch(f"{MOD}.align_reads_with_minimap2") as align,
    ):
        out = simulate_ont_fragment_pipeline(
            config, str(diploid), output_config=OutputConfig(out_dir=tmp_path, out_base="x")
        )
    return out, sim, align


def test_defaults_mix_strands_and_derive_read_count(tmp_path: Path, diploid: Path) -> None:
    config = {"read_simulation": {"coverage": 10}, "ont_fragment_params": {"length_median": 1000}}
    out, sim, align = _run(tmp_path, diploid, config)
    molecules = sim.call_args.args[0]
    assert len(molecules) == 60  # 10x * 3000 bp * 2 haplotypes / 1000 bp
    assert {m.strand for m in molecules} == {"+", "-"}
    assert out.endswith("_ont_fragments.fastq") and not align.called


def test_profile_models_and_alignment(tmp_path: Path, diploid: Path) -> None:
    profile = tmp_path / "p.json"
    profile.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "name": "g",
                "platform": "ont",
                "molecules": {"forward_frac": 1.0},
                "fragments": {"length_median": 500, "length_sigma": 0.3},
            }
        )
    )
    config = {"read_model": {"profile": str(profile)}, "ont_fragment_params": {"n_reads": 25}}
    with (
        patch(f"{MOD}.simulate_molecule_reads", return_value=0) as sim,
        patch(f"{MOD}.create_pipeline_metadata"),
        patch(f"{MOD}.align_reads_with_minimap2", return_value="x.bam") as align,
    ):
        out = simulate_ont_fragment_pipeline(
            config,
            str(diploid),
            human_reference="/ref.fa",
            output_config=OutputConfig(out_dir=tmp_path, out_base="x"),
        )
    molecules = sim.call_args.args[0]
    assert len(molecules) == 25 and all(m.strand == "+" for m in molecules)
    assert out == "x.bam" and align.call_args.kwargs["preset"] == "map-ont"

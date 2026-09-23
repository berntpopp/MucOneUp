"""End-to-end truth-tracked simulation with the built-in empirical ONT profiles.

The empirical sequencer needs no external tools, so these tests run the real
pipelines (haplotypes from the repository config, primer extraction, PCR split,
molecule model, error channel, relabelling, truth manifest) without mocks.
"""

from __future__ import annotations

import gzip
import json
from collections import Counter
from pathlib import Path

import pytest

from muc_one_up.read_simulator.ont_amplicon_pipeline import simulate_ont_amplicon_pipeline
from muc_one_up.read_simulator.ont_fragment_pipeline import simulate_ont_fragment_pipeline
from muc_one_up.read_simulator.output_config import OutputConfig
from muc_one_up.read_simulator.read_profiles import apply_read_profile, load_read_profile
from muc_one_up.simulate import simulate_from_chains
from muc_one_up.type_defs import RepeatUnit

REPO_CONFIG = Path(__file__).resolve().parents[2] / "config.json"
pytestmark = pytest.mark.integration


def _chain(n_x: int) -> list[RepeatUnit]:
    units = ["1", "2", "3", "4", "5", "C", *["X"] * n_x, "6", "7", "8", "9"]
    return [RepeatUnit.from_str(u) for u in units]


@pytest.fixture(scope="module")
def config() -> dict:
    return json.loads(REPO_CONFIG.read_text())


@pytest.fixture
def diploid_fasta(tmp_path: Path, config: dict) -> Path:
    results = simulate_from_chains([_chain(20), _chain(32)], config)
    path = tmp_path / "s.001.simulated.fa"
    path.write_text("".join(f">haplotype_{i}\n{hr.sequence}\n" for i, hr in enumerate(results, 1)))
    return path


def _truth(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"), strict=True)) for line in handle]


def _fastq_names(path: Path) -> list[str]:
    return [line[1:].strip() for i, line in enumerate(path.read_text().splitlines()) if i % 4 == 0]


def test_ont_amplicon_builtin_profile(tmp_path: Path, config: dict, diploid_fasta: Path) -> None:
    cfg = apply_read_profile(config, load_read_profile("ont_r10_sup_amplicon_v1"))
    cfg["read_simulation"] = {"coverage": 600}
    cfg["ont_amplicon_params"] = {**cfg.get("ont_amplicon_params", {}), "seed": 5}
    out = simulate_ont_amplicon_pipeline(
        cfg, str(diploid_fasta), output_config=OutputConfig(out_dir=tmp_path, out_base="amp")
    )
    rows = _truth(tmp_path / "amp_read_truth.tsv.gz")
    names = _fastq_names(Path(out))
    assert names == [r["read_id"] for r in rows] and len(set(names)) == len(names)
    kinds = Counter(r["kind"] for r in rows)
    assert kinds["smear"] > 0 and kinds["offtarget"] > 0
    strands = Counter(r["strand"] for r in rows)
    assert 0.4 < strands["+"] / sum(strands.values()) < 0.6
    assert sum(int(r["n_hp_edits"]) for r in rows) > 0


def test_ont_fragments_builtin_profile(tmp_path: Path, config: dict, diploid_fasta: Path) -> None:
    cfg = apply_read_profile(config, load_read_profile("ont_r10_genomic_v1"))
    cfg["read_simulation"] = {}
    cfg["ont_fragment_params"] = {"n_reads": 120, "length_median": 3000, "seed": 2}
    out = simulate_ont_fragment_pipeline(
        cfg, str(diploid_fasta), output_config=OutputConfig(out_dir=tmp_path, out_base="gen")
    )
    rows = _truth(tmp_path / "gen_read_truth.tsv.gz")
    assert len(rows) == 120 and _fastq_names(Path(out)) == [r["read_id"] for r in rows]
    assert {r["hap"] for r in rows} == {"1", "2"} and {r["kind"] for r in rows} == {"fragment"}

"""Shared fixtures for E2E pipeline tests.

These tests require real bioinformatics tools installed locally.
Gated by @pytest.mark.e2e and @pytest.mark.requires_tools.
"""

import json
import shutil
from pathlib import Path
from typing import Any

import pytest

from muc_one_up.simulate import simulate_from_chains
from muc_one_up.type_defs import RepeatUnit


@pytest.fixture
def e2e_config_base() -> dict[str, Any]:
    """Base config with tool paths discovered from PATH."""
    tools: dict[str, str] = {}
    for tool in [
        "samtools",
        "bwa",
        "minimap2",
        "reseq",
        "faToTwoBit",
        "pblat",
        "ccs",
        "nanosim",
    ]:
        path = shutil.which(tool)
        if path:
            tools[tool] = path
    # pbsim3 may be installed as "pbsim"
    for name in ["pbsim3", "pbsim"]:
        path = shutil.which(name)
        if path:
            tools["pbsim3"] = path
            break
    return {"tools": tools}


@pytest.fixture
def e2e_diploid_fasta(tmp_path, minimal_config) -> Path:
    """Generate a small diploid FASTA from real MUC1 repeat chains."""
    chains = [
        [RepeatUnit.from_str(s) for s in ["1", "2", "X", "B", "6", "7", "8", "9"]],
        [RepeatUnit.from_str(s) for s in ["1", "2", "A", "B", "6p", "7", "8", "9"]],
    ]
    results = simulate_from_chains(chains, minimal_config)

    fasta_path = tmp_path / "e2e_diploid.fa"
    with open(fasta_path, "w") as f:
        for hr in results:
            f.write(f">{hr.name}\n{hr.sequence}\n")
    return fasta_path


@pytest.fixture
def e2e_simulation_stats(tmp_path, e2e_diploid_fasta) -> Path:
    """Write minimal simulation_stats.json companion file."""
    stats = {
        "haplotypes": [
            {"name": "haplotype_1", "chain": "1-2-X-B-6-7-8-9"},
            {"name": "haplotype_2", "chain": "1-2-A-B-6p-7-8-9"},
        ],
        "assembly": "hg38",
    }
    stats_path = tmp_path / "e2e_diploid.simulation_stats.json"
    stats_path.write_text(json.dumps(stats))
    return stats_path

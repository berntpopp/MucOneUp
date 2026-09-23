"""CLI tests for --read-profile, amplicon --track-read-source and --no-align.

Invokes the real Click commands, patching only the batch simulator, so the
config-building logic in the commands is exercised end to end.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_up.cli.click_main import cli

PROFILE = {
    "schema_version": 1,
    "name": "cli_test_ont",
    "platform": "ont",
    "description": "cli test",
    "provenance": {"source": "test"},
    "config_overrides": {
        "ont_amplicon_params": {"accuracy_mean": 0.985, "difference_ratio": "33:31:36"},
        "amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}},
    },
    "molecules": {"forward_frac": 0.5},
}
BASE_CONFIG = {
    "tools": {"pbsim3": "pbsim", "samtools": "samtools", "minimap2": "minimap2"},
    "read_simulation": {"human_reference": "/ref/hg38.fa"},
    "amplicon_params": {"forward_primer": "ACGTACGT", "reverse_primer": "TTGGCCAA"},
    "ont_amplicon_params": {"model_file": "/m/QSHMM-ONT-HQ.model"},
    "pacbio_params": {"model_type": "errhmm", "model_file": "/m/ERRHMM-SEQUEL.model"},
}


def _invoke(tmp_path: Path, command: list[str], config: dict | None = None):
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config or BASE_CONFIG))
    fasta = tmp_path / "in.fa"
    fasta.write_text(">s\nACGT\n")
    captured: dict = {}

    def fake_batch(cfg, *args, **kwargs):
        captured["config"] = cfg

    with patch("muc_one_up.cli.commands.reads._run_batch_simulation", side_effect=fake_batch):
        result = CliRunner().invoke(
            cli, ["--config", str(config_path), "reads", *command, str(fasta)]
        )
    return result, captured.get("config")


def _profile(tmp_path: Path, **changes) -> str:
    path = tmp_path / "profile.json"
    path.write_text(json.dumps({**PROFILE, **changes}))
    return str(path)


def test_amplicon_track_read_source_is_accepted(tmp_path):
    result, config = _invoke(tmp_path, ["amplicon", "--platform", "ont", "--track-read-source"])
    assert result.exit_code == 0, result.output
    assert config["read_simulation"]["track_read_source"] is True


def test_read_profile_overlays_config(tmp_path):
    result, config = _invoke(
        tmp_path, ["amplicon", "--platform", "ont", "--read-profile", _profile(tmp_path)]
    )
    assert result.exit_code == 0, result.output
    assert config["ont_amplicon_params"]["difference_ratio"] == "33:31:36"
    assert config["ont_amplicon_params"]["model_file"] == "/m/QSHMM-ONT-HQ.model"
    assert config["amplicon_params"]["pcr_bias"]["preset"] == "madritsch2025_r10"
    assert config["read_model"]["name"] == "cli_test_ont"


def test_cli_flags_override_profile(tmp_path):
    result, config = _invoke(
        tmp_path,
        [
            "amplicon",
            "--platform",
            "ont",
            "--read-profile",
            _profile(tmp_path),
            "--pcr-preset",
            "no_bias",
        ],
    )
    assert result.exit_code == 0, result.output
    assert config["amplicon_params"]["pcr_bias"]["preset"] == "no_bias"


def test_profile_platform_mismatch_is_rejected(tmp_path):
    result, _ = _invoke(tmp_path, ["amplicon", "--read-profile", _profile(tmp_path)])
    assert result.exit_code != 0
    assert "platform" in result.output


def test_profile_named_in_config_is_applied(tmp_path):
    config = {**BASE_CONFIG, "read_model": {"profile": _profile(tmp_path)}}
    result, merged = _invoke(tmp_path, ["amplicon", "--platform", "ont"], config)
    assert result.exit_code == 0, result.output
    assert merged["ont_amplicon_params"]["accuracy_mean"] == 0.985


def test_no_align_drops_reference(tmp_path):
    for command in (["amplicon", "--platform", "ont"], ["ont"], ["pacbio"]):
        result, config = _invoke(tmp_path, [*command, "--no-align"])
        assert result.exit_code == 0, result.output
        assert "human_reference" not in config["read_simulation"]
        assert config["read_simulation"]["skip_alignment"] is True  # honoured by NanoSim (#116)


def test_reference_kept_without_no_align(tmp_path):
    _result, config = _invoke(tmp_path, ["amplicon", "--platform", "ont"])
    assert config["read_simulation"]["human_reference"] == "/ref/hg38.fa"


def test_profiles_command_lists_builtins():
    result = CliRunner().invoke(cli, ["reads", "profiles"])
    assert result.exit_code == 0, result.output


def test_ont_pbsim3_fragments_configures_simulator(tmp_path):
    result, config = _invoke(
        tmp_path,
        [
            "ont",
            "--simulator",
            "pbsim3-fragments",
            "--read-profile",
            _profile(tmp_path),
            "--n-reads",
            "300",
            "--read-length-median",
            "6000",
            "--seed",
            "4",
        ],
    )
    assert result.exit_code == 0, result.output
    assert config["read_simulation"]["simulator"] == "ont-fragments"
    assert config["ont_fragment_params"] == {"n_reads": 300, "length_median": 6000.0, "seed": 4}
    assert config["read_model"]["name"] == "cli_test_ont"


def test_ont_read_profile_requires_fragment_simulator(tmp_path):
    result, _ = _invoke(tmp_path, ["ont", "--read-profile", _profile(tmp_path)])
    assert result.exit_code != 0 and "pbsim3-fragments" in result.output

"""Tests for ONT CLI config key mapping.

These tests invoke the actual ont() Click command through CliRunner,
patching only the config loader and batch simulator. This ensures the
real config-building logic in ont() is exercised.
"""

import json
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


def _invoke_ont(tmp_path, config_dict, cli_args=None):
    """Invoke 'muconeup reads ont' and return the config passed to the simulator.

    Writes config_dict to a real JSON file (Click validates --config exists),
    patches load_config_raw to return config_dict, and captures
    the config that _run_batch_simulation receives.
    """
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config_dict))

    # Also need a fake input FASTA
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">seq\nACGT\n")

    captured = {}

    def fake_load_config_raw(path):
        return json.loads(config_path.read_text())

    def fake_batch_sim(config, *args, **kwargs):
        captured["config"] = config

    runner = CliRunner()
    with (
        patch(
            "muc_one_up.config.load_config_raw",
            side_effect=fake_load_config_raw,
        ),
        patch(
            "muc_one_up.cli.commands.reads._run_batch_simulation",
            side_effect=fake_batch_sim,
        ),
    ):
        args = ["--config", str(config_path), "reads", "ont"]
        if cli_args:
            args.extend(cli_args)
        args.append(str(fasta_path))
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, f"CLI failed: {result.output}"
    return captured["config"]


def test_ont_cli_coverage_lands_in_nanosim_params(tmp_path):
    """CLI --coverage must propagate to nanosim_params.coverage."""
    base = {"nanosim_params": {"training_data_path": "/fake"}}
    config = _invoke_ont(tmp_path, base, ["--coverage", "75"])
    assert config["nanosim_params"]["coverage"] == 75


def test_ont_config_file_coverage_not_clobbered(tmp_path):
    """Config-file nanosim_params.coverage must survive when CLI omits --coverage."""
    base = {
        "nanosim_params": {"training_data_path": "/fake", "coverage": 200},
        "read_simulation": {"coverage": 30},
    }
    config = _invoke_ont(tmp_path, base)
    # Config-file value of 200 must NOT be clobbered by the default 30
    assert config["nanosim_params"]["coverage"] == 200


def test_ont_min_read_length_uses_correct_key(tmp_path):
    """CLI --min-read-length must map to nanosim_params.min_read_length."""
    base = {"nanosim_params": {"training_data_path": "/fake", "coverage": 50}}
    config = _invoke_ont(tmp_path, base, ["--min-read-length", "500"])
    assert config["nanosim_params"]["min_read_length"] == 500
    assert "min_len" not in config["nanosim_params"]


def test_ont_default_coverage_propagated_when_missing(tmp_path):
    """When neither CLI nor config provides ONT coverage, default propagates."""
    base = {"nanosim_params": {"training_data_path": "/fake"}}
    config = _invoke_ont(tmp_path, base)
    # _setup_read_config defaults to 30 in read_simulation; should propagate
    assert config["nanosim_params"]["coverage"] == 30

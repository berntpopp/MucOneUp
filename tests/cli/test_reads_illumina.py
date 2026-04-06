"""Tests for reads illumina --read-number CLI option."""

import json
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


def _invoke_illumina(tmp_path, read_number=None, coverage=None):
    """Invoke reads illumina and return the config passed to the simulator."""
    config_dict = {
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref.fa",
            "reseq_model": "/model.reseq",
            "read_number": 50000,
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config_dict))

    fasta = tmp_path / "input.fa"
    fasta.write_text(">seq\nACGT\n")

    captured = {}

    def fake_load(path):
        return json.loads(config_path.read_text())

    def fake_batch(config, *args, **kwargs):
        captured["config"] = config

    runner = CliRunner()
    with (
        patch("muc_one_up.config.load_config_raw", side_effect=fake_load),
        patch(
            "muc_one_up.cli.commands.reads._run_batch_simulation",
            side_effect=fake_batch,
        ),
    ):
        args = ["--config", str(config_path), "reads", "illumina"]
        if read_number is not None:
            args.extend(["--read-number", str(read_number)])
        if coverage is not None:
            args.extend(["--coverage", str(coverage)])
        args.append(str(fasta))
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, f"CLI failed: {result.output}"
    return captured["config"]


def test_default_read_number_from_config(tmp_path):
    """Without --read-number, config value is preserved."""
    config = _invoke_illumina(tmp_path)
    assert config["read_simulation"]["read_number"] == 50000


def test_read_number_overrides_config(tmp_path):
    """--read-number overrides the config value."""
    config = _invoke_illumina(tmp_path, read_number=200000)
    assert config["read_simulation"]["read_number"] == 200000


def test_read_number_with_coverage(tmp_path):
    """--read-number and --coverage can be used together."""
    config = _invoke_illumina(tmp_path, read_number=500000, coverage=1000)
    assert config["read_simulation"]["read_number"] == 500000
    assert config["read_simulation"]["coverage"] == 1000

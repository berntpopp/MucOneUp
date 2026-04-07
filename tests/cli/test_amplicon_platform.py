"""Tests for reads amplicon --platform routing."""

import json
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


def _invoke_amplicon(tmp_path, platform=None, extra_args=None):
    """Invoke reads amplicon and return the config passed to the simulator."""
    model_file = tmp_path / "test.model"
    model_file.write_text("model")

    config_dict = {
        "amplicon_params": {
            "forward_primer": "ACGT",
            "reverse_primer": "TGCA",
        },
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": str(model_file),
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
        args = ["--config", str(config_path), "reads", "amplicon"]
        if platform:
            args.extend(["--platform", platform])
        if extra_args:
            args.extend(extra_args)
        args.append(str(fasta))
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, f"CLI failed: {result.output}"
    return captured["config"]


def test_default_platform_is_pacbio(tmp_path):
    """Without --platform, simulator should be 'amplicon' (PacBio)."""
    config = _invoke_amplicon(tmp_path)
    assert config["read_simulation"]["simulator"] == "amplicon"


def test_platform_ont_sets_ont_amplicon(tmp_path):
    """--platform ont should set simulator to 'ont-amplicon'."""
    config = _invoke_amplicon(tmp_path, platform="ont")
    assert config["read_simulation"]["simulator"] == "ont-amplicon"


def test_ont_platform_uses_ont_amplicon_params(tmp_path):
    """--platform ont should populate ont_amplicon_params, not pacbio_params."""
    config = _invoke_amplicon(tmp_path, platform="ont")
    assert "ont_amplicon_params" in config
    # Default ONT model should be set
    assert "ONT" in config["ont_amplicon_params"]["model_file"].upper()


def test_assay_type_set_for_both_platforms(tmp_path):
    """Both platforms should set assay_type='amplicon'."""
    for platform in [None, "ont"]:
        config = _invoke_amplicon(tmp_path, platform=platform)
        assert config["read_simulation"]["assay_type"] == "amplicon"

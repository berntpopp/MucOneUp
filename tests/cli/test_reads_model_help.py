"""CLI help must only reference pbsim3 models that pbsim3 actually ships (#110)."""

import pytest
from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


@pytest.mark.parametrize("command", ["pacbio", "amplicon"])
def test_help_does_not_reference_missing_qshmm_sequel_model(command):
    result = CliRunner().invoke(cli, ["reads", command, "--help"])

    assert result.exit_code == 0
    assert "QSHMM-SEQUEL" not in result.output
    assert "--model-type errhmm --model-file /models/ERRHMM-SEQUEL.model" in result.output

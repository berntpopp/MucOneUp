"""Tests for PacBio parameter consolidation in CLI."""

import logging

import click
import pytest


class TestApplyPacbioParams:
    """Test the shared PacBio parameter helper."""

    def test_initializes_pacbio_params(self):
        """Creates pacbio_params section if missing."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        _apply_pacbio_params(config, model_type="qshmm", model_file="/model", seed=None)
        assert "pacbio_params" in config
        assert config["pacbio_params"]["model_type"] == "qshmm"
        assert config["pacbio_params"]["model_file"] == "/model"

    def test_only_overrides_non_none(self):
        """Only sets params that are not None."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config = {"pacbio_params": {"model_type": "errhmm", "threads": 8}}
        _apply_pacbio_params(
            config, model_type=None, model_file="/model", seed=None, threads=None
        )
        assert config["pacbio_params"]["model_type"] == "errhmm"  # not overwritten
        assert config["pacbio_params"]["model_file"] == "/model"
        assert config["pacbio_params"]["threads"] == 8  # not overwritten

    def test_raises_on_missing_required(self):
        """Raises ClickException when required params missing."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        with pytest.raises(click.ClickException, match="model_type"):
            _apply_pacbio_params(config, model_type=None, model_file=None, seed=None)

    def test_logs_seed(self, caplog):
        """Logs seed when provided."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        with caplog.at_level(logging.INFO):
            _apply_pacbio_params(
                config, model_type="qshmm", model_file="/m", seed=42
            )
        assert "42" in caplog.text
        assert config["pacbio_params"]["seed"] == 42

    def test_sets_optional_params(self):
        """Sets optional params when provided."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        _apply_pacbio_params(
            config, model_type="qshmm", model_file="/m", seed=None,
            threads=8, pass_num=5, min_passes=3, min_rq=0.999,
        )
        assert config["pacbio_params"]["threads"] == 8
        assert config["pacbio_params"]["pass_num"] == 5
        assert config["pacbio_params"]["min_passes"] == 3
        assert config["pacbio_params"]["min_rq"] == 0.999

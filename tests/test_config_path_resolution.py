"""Tests for config path resolution (#80)."""

import json

from muc_one_up.config import load_config_raw


def test_resolves_relative_paths_from_config_dir(tmp_path, monkeypatch):
    """Relative paths in config are resolved relative to config file directory."""
    ref_dir = tmp_path / "reference"
    ref_dir.mkdir()
    model = ref_dir / "test.model"
    model.write_text("model")

    config = {
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": "reference/test.model",
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    # Load from a different working directory
    monkeypatch.chdir("/tmp")
    loaded = load_config_raw(str(config_path))
    assert loaded["pacbio_params"]["model_file"] == str(model)


def test_preserves_absolute_paths(tmp_path):
    """Absolute paths in config are not modified."""
    config = {
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": "/absolute/path/to/model.file",
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    loaded = load_config_raw(str(config_path))
    assert loaded["pacbio_params"]["model_file"] == "/absolute/path/to/model.file"


def test_preserves_paths_that_exist_from_cwd(tmp_path, monkeypatch):
    """Paths that already resolve from cwd are not changed."""
    config = {
        "read_simulation": {
            "human_reference": "config.json",
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    monkeypatch.chdir(tmp_path)
    loaded = load_config_raw(str(config_path))
    assert loaded["read_simulation"]["human_reference"] == "config.json"

import os
import pytest
import json
from muc_one_up.config import load_config

def test_load_config_valid(tmp_path):
    """Test that a valid config file loads without error."""
    config_data = {
        "repeats": {"1": "ABC"},
        "constants": {"left": "AAA", "right": "TTT"},
        "probabilities": {"1": {"2": 1.0}},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 5,
            "max_repeats": 10,
            "mean_repeats": 7
        },
        "mutations": {}
    }
    config_file = tmp_path / "config.json"
    with open(config_file, "w") as fh:
        json.dump(config_data, fh)

    config = load_config(str(config_file))
    assert "repeats" in config
    assert config["repeats"]["1"] == "ABC"
    assert config["constants"]["left"] == "AAA"

def test_load_config_file_not_found():
    """Test that loading a non-existent file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_config("no_such_file.json")

def test_load_config_bad_json(tmp_path):
    """Test that a malformed JSON file raises JSONDecodeError."""
    bad_file = tmp_path / "bad.json"
    with open(bad_file, "w") as fh:
        fh.write("{ invalid json }")

    with pytest.raises(json.JSONDecodeError):
        load_config(str(bad_file))

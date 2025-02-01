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
            "mean_repeats": 7,
            "median_repeats": 7
        },
        "mutations": {},
        "tools": {
            "reseq": "dummy_reseq",
            "faToTwoBit": "dummy_faToTwoBit",
            "samtools": "dummy_samtools",
            "pblat": "dummy_pblat",
            "bwa": "dummy_bwa"
        },
        "read_simulation": {
            "reseq_model": "dummy_reseq_model",
            "sample_bam": "dummy_sample_bam",
            "human_reference": "dummy_human_reference",
            "read_number": 100,
            "fragment_size": 150,
            "fragment_sd": 20,
            "min_fragment": 30,
            "threads": 1
        }
    }
    config_file = tmp_path / "config.json"
    with open(config_file, "w") as fh:
        json.dump(config_data, fh)

    config = load_config(str(config_file))
    assert "repeats" in config
    assert config["repeats"]["1"] == "ABC"
    assert config["constants"]["left"] == "AAA"

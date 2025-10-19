import json
from pathlib import Path

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
            "median_repeats": 7,
        },
        "mutations": {},
        "tools": {
            "reseq": "dummy_reseq",
            "faToTwoBit": "dummy_faToTwoBit",
            "samtools": "dummy_samtools",
            "pblat": "dummy_pblat",
            "bwa": "dummy_bwa",
        },
        "read_simulation": {
            "reseq_model": "dummy_reseq_model",
            "sample_bam": "dummy_sample_bam",
            "human_reference": "dummy_human_reference",
            "read_number": 100,
            "fragment_size": 150,
            "fragment_sd": 20,
            "min_fragment": 30,
            "threads": 1,
        },
    }
    config_file = tmp_path / "config.json"
    with Path(config_file).open("w") as fh:
        json.dump(config_data, fh)

    config = load_config(str(config_file))
    assert "repeats" in config
    assert config["repeats"]["1"] == "ABC"
    # After normalization, flat constants are converted to nested format (default: hg38)
    ref_assembly = config.get("reference_assembly", "hg38")
    assert config["constants"][ref_assembly]["left"] == "AAA"


def test_read_simulation_accepts_seed(tmp_path):
    """Test that read_simulation section accepts seed parameter."""
    config_data = {
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
            "seed": 42,  # Should be valid
        },
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    # Should not raise ValidationError
    config = load_config(str(config_file))
    assert config["read_simulation"]["seed"] == 42


def test_nanosim_params_accepts_seed(tmp_path):
    """Test that nanosim_params section accepts seed parameter."""
    config_data = {
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
        },
        "nanosim_params": {
            "training_data_path": "/path/to/model",
            "coverage": 30,
            "seed": 12345,  # Should be valid
        },
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    # Should not raise ValidationError
    config = load_config(str(config_file))
    assert config["nanosim_params"]["seed"] == 12345


def test_seed_can_be_null(tmp_path):
    """Test that seed can be explicitly null (for JSON compatibility)."""
    config_data = {
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
            "seed": None,  # Explicit null
        },
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    # Should not raise ValidationError
    config = load_config(str(config_file))
    assert config["read_simulation"]["seed"] is None

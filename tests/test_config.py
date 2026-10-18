import json
from pathlib import Path

import pytest
from jsonschema import ValidationError

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


def test_load_config_accepts_root_seed(tmp_path):
    """Test that root-level seed is accepted by load_config."""
    config_data = {
        "seed": 42,
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
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    loaded = load_config(str(config_file))
    assert loaded["seed"] == 42

    # Also test seed: null
    config_data["seed"] = None
    config_file.write_text(json.dumps(config_data))
    loaded = load_config(str(config_file))
    assert loaded["seed"] is None


def test_load_config_accepts_root_reference_assembly(tmp_path):
    """Test that root-level reference_assembly is accepted and validated by load_config."""
    config_data = {
        "reference_assembly": "hg38",
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
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    loaded = load_config(str(config_file))
    assert loaded["reference_assembly"] == "hg38"

    config_data["reference_assembly"] = "hg19"
    config_file.write_text(json.dumps(config_data))
    loaded = load_config(str(config_file))
    assert loaded["reference_assembly"] == "hg19"

    # Invalid assembly should raise ValidationError
    config_data["reference_assembly"] = "hg37"
    config_file.write_text(json.dumps(config_data))
    with pytest.raises(ValidationError):
        load_config(str(config_file))


def test_load_config_resolves_relative_paths(tmp_path, monkeypatch):
    """Test that load_config resolves relative paths against config directory."""
    ref_file = tmp_path / "reference.fa"
    ref_file.write_text(">ref\nACGT\n")

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
            "human_reference": "reference.fa",
            "threads": 4,
        },
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    monkeypatch.chdir("/tmp")
    loaded = load_config(str(config_file))
    assert loaded["read_simulation"]["human_reference"] == str(ref_file.resolve())


def test_config_schema_accepts_snapshot_primary_primer_and_mutant_extension():
    """Verify CONFIG_SCHEMA accepts snapshot.primary_primer and validation.mutant_extension."""
    from jsonschema import validate

    from muc_one_up.config import CONFIG_SCHEMA

    snapshot_props = CONFIG_SCHEMA["properties"]["snapshot_validation"]["additionalProperties"][
        "properties"
    ]["snapshot"]["properties"]
    assert "primary_primer" in snapshot_props
    assert snapshot_props["primary_primer"] == {"type": "string"}

    val_props = CONFIG_SCHEMA["properties"]["snapshot_validation"]["additionalProperties"][
        "properties"
    ]["validation"]["properties"]
    assert "mutant_extension" in val_props
    assert val_props["mutant_extension"] == {"type": "object"}

    sample_config = {
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
        "snapshot_validation": {
            "test_assay": {
                "pcr": {
                    "forward_primer": "ACGT",
                    "reverse_primer": "TGCA",
                    "reverse_needs_rc": False,
                    "max_products": 1,
                    "size_range": {"min": 100, "max": 500},
                },
                "digest": {
                    "enzyme": "EcoRI",
                    "recognition_site": "GAATTC",
                },
                "snapshot": {
                    "primary_primer": "ACTG",
                    "primers": {"p1": "ACTG"},
                    "fluorophore_map": {
                        "dC": {"color": "blue", "dye": "FAM"},
                    },
                },
                "validation": {
                    "mutant_pattern": ".*",
                    "normal_pattern": ".*",
                    "mutant_extension": {"bases": "C"},
                },
            }
        },
    }
    validate(instance=sample_config, schema=CONFIG_SCHEMA)


def test_load_config_raw_rejects_non_dict_json(tmp_path):
    """load_config_raw should raise ConfigurationError when JSON is not an object."""
    from muc_one_up.config import load_config_raw
    from muc_one_up.exceptions import ConfigurationError

    array_config = tmp_path / "array.json"
    array_config.write_text("[1, 2, 3]")

    with pytest.raises(ConfigurationError, match="must contain a JSON object"):
        load_config_raw(str(array_config))


def test_validate_simulation_config_valid(minimal_config):
    """validate_simulation_config succeeds on valid minimal_config."""
    from muc_one_up.config import validate_simulation_config

    validate_simulation_config(minimal_config)


def test_validate_simulation_config_missing_required_section():
    """validate_simulation_config raises ConfigurationError when a required section is missing."""
    from muc_one_up.config import validate_simulation_config
    from muc_one_up.exceptions import ConfigurationError

    cfg = {"repeats": {"1": "ATCG"}, "probabilities": {}}
    with pytest.raises(ConfigurationError, match="Missing or null required configuration section"):
        validate_simulation_config(cfg)


def test_validate_simulation_config_invalid_schema():
    """validate_simulation_config raises ConfigurationError when section schema is invalid."""
    from muc_one_up.config import validate_simulation_config
    from muc_one_up.exceptions import ConfigurationError

    cfg = {
        "repeats": "not_a_dict",  # invalid
        "constants": {"hg38": {"left": "A", "right": "T"}},
        "probabilities": {},
    }
    with pytest.raises(ConfigurationError, match="Invalid configuration in section 'repeats'"):
        validate_simulation_config(cfg)

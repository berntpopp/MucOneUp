# muc_one_up/config.py

import json
import os
import logging
from typing import Any, Dict

from jsonschema import validate, ValidationError

# Define a JSON schema that matches the current config structure.
CONFIG_SCHEMA: Dict[str, Any] = {
    "type": "object",
    "properties": {
        "repeats": {"type": "object", "additionalProperties": {"type": "string"}},
        "nanosim_params": {
            "type": "object",
            "properties": {
                "training_data_path": {"type": "string"},
                "coverage": {"type": "number"},
                "num_threads": {"type": "number"},
                "min_read_length": {"type": "number"},
                "max_read_length": {"type": "number"},
                "other_options": {"type": "string"},
            },
            "required": ["training_data_path", "coverage"],
            "additionalProperties": False,
        },
        "constants": {
            "type": "object",
            "properties": {"left": {"type": "string"}, "right": {"type": "string"}},
            "required": ["left", "right"],
            "additionalProperties": False,
        },
        "probabilities": {
            "type": "object",
            "additionalProperties": {
                "type": "object",
                "additionalProperties": {"type": "number"},
            },
        },
        "length_model": {
            "type": "object",
            "properties": {
                "distribution": {"type": "string"},
                "min_repeats": {"type": "number"},
                "max_repeats": {"type": "number"},
                "mean_repeats": {"type": "number"},
                "median_repeats": {"type": "number"},
            },
            "required": [
                "distribution",
                "min_repeats",
                "max_repeats",
                "mean_repeats",
                "median_repeats",
            ],
            "additionalProperties": False,
        },
        "mutations": {
            "type": "object",
            "additionalProperties": {
                "type": "object",
                "properties": {
                    "allowed_repeats": {"type": "array", "items": {"type": "string"}},
                    "strict_mode": {"type": "boolean"},
                    "changes": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "type": {
                                    "type": "string",
                                    "enum": [
                                        "insert",
                                        "delete",
                                        "replace",
                                        "delete_insert",
                                    ],
                                },
                                "start": {"type": "number"},
                                "end": {"type": "number"},
                                "sequence": {"type": "string"},
                            },
                            "required": ["type", "start", "end"],
                            "additionalProperties": False,
                        },
                    },
                },
                "required": ["allowed_repeats", "changes"],
                "additionalProperties": False,
            },
        },
        "tools": {
            "type": "object",
            "properties": {
                "reseq": {"type": "string"},
                "faToTwoBit": {"type": "string"},
                "samtools": {"type": "string"},
                "pblat": {"type": "string"},
                "bwa": {"type": "string"},
                "nanosim": {"type": "string"},
                "minimap2": {"type": "string"},
            },
            "required": ["samtools"],
            "additionalProperties": False,
        },
        "read_simulation": {
            "type": "object",
            "properties": {
                "simulator": {"type": "string", "enum": ["illumina", "ont"]},
                "reseq_model": {"type": "string"},
                "sample_bam": {"type": "string"},
                "human_reference": {"type": "string"},
                "read_number": {"type": "number"},
                "fragment_size": {"type": "number"},
                "fragment_sd": {"type": "number"},
                "min_fragment": {"type": "number"},
                "threads": {"type": "number"},
                "downsample_coverage": {"type": "number"},
                "downsample_seed": {"type": "number"},
                "reference_assembly": {"type": "string"},
                "vntr_region_hg19": {"type": "string"},
                "vntr_region_hg38": {"type": "string"},
                "aligner": {"type": "string", "enum": ["bwa", "minimap2"]},
            },
            "required": [
                "simulator",
                "human_reference",
                "threads",
            ],
            "additionalProperties": False,
        },
    },
    "required": [
        "repeats",
        "constants",
        "probabilities",
        "length_model",
        "mutations",
        "tools",
        "read_simulation",
    ],
    "additionalProperties": False,
}


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load and validate the JSON config for MucOneUp.

    :param config_path: Path to the JSON config file.
    :return: Python dict with validated configuration data.
    :raises FileNotFoundError: if the config file does not exist.
    :raises json.JSONDecodeError: if the config file is not valid JSON.
    :raises ValidationError: if the config does not conform to the required schema.
    """
    if not os.path.exists(config_path):
        logging.error("Config file not found: %s", config_path)
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as fh:
        config = json.load(fh)

    logging.debug("Configuration loaded from %s", config_path)

    try:
        validate(instance=config, schema=CONFIG_SCHEMA)

        # Extra validation for allowed_repeats in mutations
        if "repeats" in config and "mutations" in config:
            repeat_keys = set(config["repeats"].keys())
            for mutation_name, mutation_def in config["mutations"].items():
                allowed_repeats = set(mutation_def.get("allowed_repeats", []))
                invalid_repeats = allowed_repeats - repeat_keys

                if invalid_repeats:
                    msg = (
                        f"Mutation '{mutation_name}' has invalid repeats "
                        f"in allowed_repeats: {', '.join(invalid_repeats)}. "
                        f"Valid repeats are: {', '.join(repeat_keys)}"
                    )
                    logging.error(msg)
                    raise ValidationError(msg)
    except ValidationError as e:
        logging.error("Configuration validation error: %s", e.message)
        raise

    return config

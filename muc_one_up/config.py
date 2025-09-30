# muc_one_up/config.py

import json
import logging
from pathlib import Path
from typing import Any

from jsonschema import ValidationError, validate

# Define a JSON schema that matches the current config structure.
CONFIG_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "repeats": {"type": "object", "additionalProperties": {"type": "string"}},
        "nanosim_params": {
            "type": "object",
            "properties": {
                "training_data_path": {"type": "string"},
                "coverage": {"type": "number"},
                "num_threads": {"type": ["number", "null"]},
                "min_read_length": {"type": ["number", "null"]},
                "max_read_length": {"type": ["number", "null"]},
                "other_options": {"type": ["string", "null"]},
            },
            "required": ["training_data_path", "coverage"],
            "additionalProperties": False,
        },
        "constants": {
            "type": "object",
            # Allow either flat format (for backward compatibility with tests)
            # or nested hg19/hg38 format (for production)
            "additionalProperties": {
                "oneOf": [
                    {"type": "string"},  # Flat format: left, right keys
                    {  # Nested format: hg19, hg38 keys
                        "type": "object",
                        "properties": {
                            "left": {"type": "string"},
                            "right": {"type": "string"},
                            "vntr_region": {"type": "string"},
                        },
                    },
                ]
            },
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
                "sample_bam_hg19": {"type": "string"},
                "sample_bam_hg38": {"type": "string"},
                "human_reference": {"type": "string"},
                "read_number": {"type": "number"},
                "fragment_size": {"type": "number"},
                "fragment_sd": {"type": "number"},
                "min_fragment": {"type": "number"},
                "threads": {"type": "number"},
                "downsample_coverage": {"type": "number"},
                "downsample_seed": {"type": "number"},
                "downsample_mode": {"type": "string"},
                "sample_target_bed": {"type": "string"},
                "reference_assembly": {"type": "string"},
                "vntr_region_hg19": {"type": "string"},
                "vntr_region_hg38": {"type": "string"},
                "aligner": {"type": "string", "enum": ["bwa", "minimap2"]},
            },
            "required": [
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


def load_config(config_path: str) -> dict[str, Any]:
    """
    Load and validate the JSON config for MucOneUp.

    :param config_path: Path to the JSON config file.
    :return: Python dict with validated configuration data.
    :raises FileNotFoundError: if the config file does not exist.
    :raises json.JSONDecodeError: if the config file is not valid JSON.
    :raises ValidationError: if the config does not conform to the required schema.
    """
    if not Path(config_path).exists():
        logging.error("Config file not found: %s", config_path)
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with Path(config_path).open() as fh:
        config = json.load(fh)

    logging.debug("Configuration loaded from %s", config_path)

    try:
        validate(instance=config, schema=CONFIG_SCHEMA)

        # Normalize constants format: convert flat format to nested format
        if (
            "constants" in config
            and "left" in config["constants"]
            and "right" in config["constants"]
        ):
            # Check if this is the flat format (has 'left' and 'right' as top-level keys)
            # Get reference assembly, default to hg38 if not specified
            ref_assembly = config.get("reference_assembly", "hg38")
            # Convert to nested format
            old_constants = config["constants"].copy()
            config["constants"] = {ref_assembly: old_constants}
            logging.debug("Normalized flat constants format to nested format")

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

# muc_one_up/config.py
"""Configuration loading and validation for MucOneUp.

This module handles JSON configuration file parsing, schema validation, and
format normalization. Validates all required sections (repeats, constants,
probabilities, mutations, tools, read_simulation) and enforces structural
constraints via CONFIG_SCHEMA.

Key Functions:
    load_config: Load, validate, and normalize configuration from JSON file

Key Constants:
    CONFIG_SCHEMA: JSON Schema defining required configuration structure

Configuration Sections:
    - repeats: Mapping of repeat symbols (1, 2, X, A, B...) to DNA sequences
    - constants: Left/right flanking sequences for hg19 and hg38 assemblies
    - probabilities: State transition matrix for repeat chain generation
    - length_model: Distribution parameters (normal/uniform) for VNTR length
    - mutations: Named mutation definitions with allowed_repeats and changes
    - tools: Command paths for external tools (reseq, bwa, samtools...)
    - read_simulation: Parameters for Illumina pipeline (coverage, threads...)
    - nanosim_params: Parameters for ONT pipeline (model path, read lengths, split-simulation...)
    - pacbio_params: Parameters for PacBio HiFi pipeline (model type/file, coverage, pass numbers...)

Example:
    Load and access configuration::

        from muc_one_up.config import load_config

        config = load_config("config.json")
        repeats = config["repeats"]  # Dict[str, str]
        probs = config["probabilities"]  # Dict[str, Dict[str, float]]

Notes:
    - Automatically normalizes flat constants format to nested format
    - Validates mutation allowed_repeats against valid repeat symbols
    - Supports both hg19 and hg38 reference assemblies

Raises:
    FileNotFoundError: If config file doesn't exist
    json.JSONDecodeError: If config file contains invalid JSON
    ValidationError: If config doesn't conform to CONFIG_SCHEMA
"""

import json
import logging
from pathlib import Path
from typing import Any

from jsonschema import ValidationError, validate

#: JSON Schema for validating MucOneUp configuration files.
#:
#: Defines the complete structure and validation rules for configuration
#: dictionaries. All sections are required except where specified.
#:
#: Required Sections:
#:     repeats: Mapping of repeat symbols (str) to DNA sequences (str)
#:         Example: {"1": "ACGACT", "2": "CAGACT", "X": "TCGACT", ...}
#:
#:     constants: Flanking sequences for each reference assembly
#:         Nested format: {assembly: {left: str, right: str, vntr_region: str}}
#:         Supported assemblies: hg19, hg38
#:
#:     probabilities: State transition matrix for repeat chain generation
#:         Format: {current_symbol: {next_symbol: probability, ...}, ...}
#:         Example: {"1": {"2": 0.7, "7": 0.3}, "2": {"7": 1.0}}
#:
#:     length_model: Distribution parameters for VNTR length sampling
#:         Required keys: distribution, min_repeats, max_repeats,
#:                       mean_repeats, median_repeats
#:
#:     mutations: Named mutation definitions
#:         Each mutation has:
#:         - allowed_repeats: List of valid repeat symbols for this mutation
#:         - strict_mode: Boolean (optional, default: false)
#:         - changes: List of operations with type, start, end, sequence
#:         Operation types: insert, delete, replace, delete_insert
#:
#:     tools: Command paths for external executables
#:         Required: samtools
#:         Optional: reseq, faToTwoBit, pblat, bwa, nanosim, minimap2, pbsim3, ccs
#:
#:     read_simulation: Parameters for Illumina read simulation
#:         Required: human_reference, threads
#:         Simulator type: illumina or ont
#:
#:     nanosim_params: Parameters for Oxford Nanopore simulation (optional)
#:         Required: training_data_path, coverage
#:         Optional: num_threads, min_read_length, max_read_length, other_options, seed,
#:                   correction_factor (default: 0.325),
#:                   enable_split_simulation (default: True),
#:                   enable_coverage_correction (default: True)
#:
#:     pacbio_params: Parameters for PacBio HiFi simulation (optional)
#:         Required: model_type, model_file, coverage, pass_num, min_passes, min_rq
#:         Optional: pbsim3_cmd, ccs_cmd, threads, seed,
#:                   accuracy_mean, accuracy_sd, accuracy_min,
#:                   length_mean, length_sd, length_min, length_max
#:         Enums: model_type must be "qshmm" or "errhmm"
#:         Ranges: coverage (0.1-10000), pass_num (2-50), min_passes (1-50),
#:                min_rq (0.0-1.0), accuracy_* (0.0-1.0)
#:
#: Example:
#:     Accessing schema in code::
#:
#:         from muc_one_up.config import CONFIG_SCHEMA
#:         from jsonschema import validate
#:
#:         validate(instance=config_dict, schema=CONFIG_SCHEMA)
#:
#: Notes:
#:     - Constants support both flat format (backward compatibility) and
#:       nested format (hg19/hg38)
#:     - Mutations are validated to ensure allowed_repeats are valid symbols
#:     - All numeric fields are validated for type and range
CONFIG_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "repeats": {"type": "object", "additionalProperties": {"type": "string"}},
        "reference_genomes": {
            "type": "object",
            "patternProperties": {
                "^[a-zA-Z0-9_]+$": {  # Assembly name pattern
                    "type": "object",
                    "properties": {
                        "fasta_path": {"type": "string"},
                        "vntr_region": {"type": "string"},
                        "display_name": {"type": "string"},
                        "source_url": {"type": "string"},
                    },
                    "required": ["fasta_path", "vntr_region"],
                    "additionalProperties": False,
                }
            },
            "additionalProperties": False,
        },
        "nanosim_params": {
            "type": "object",
            "properties": {
                "training_data_path": {"type": "string"},
                "coverage": {"type": "number"},
                "num_threads": {"type": ["number", "null"]},
                "min_read_length": {"type": ["number", "null"]},
                "max_read_length": {"type": ["number", "null"]},
                "other_options": {"type": ["string", "null"]},
                "seed": {"type": ["number", "null"]},
                "correction_factor": {"type": ["number", "null"]},
                "enable_split_simulation": {"type": ["boolean", "null"]},
                "enable_coverage_correction": {"type": ["boolean", "null"]},
            },
            "required": ["training_data_path", "coverage"],
            "additionalProperties": False,
        },
        "pacbio_params": {
            "type": "object",
            "properties": {
                "pbsim3_cmd": {"type": "string"},
                "ccs_cmd": {"type": "string"},
                "model_type": {"type": "string", "enum": ["qshmm", "errhmm"]},
                "model_file": {"type": "string"},
                "coverage": {"type": "number", "minimum": 0.1, "maximum": 10000},
                "pass_num": {"type": "number", "minimum": 2, "maximum": 50},
                "min_passes": {"type": "number", "minimum": 1, "maximum": 50},
                "min_rq": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                "threads": {"type": "number", "minimum": 1},
                "seed": {"type": ["number", "null"]},
                "accuracy_mean": {"type": ["number", "null"], "minimum": 0.0, "maximum": 1.0},
                "accuracy_sd": {"type": ["number", "null"], "minimum": 0.0, "maximum": 1.0},
                "accuracy_min": {"type": ["number", "null"], "minimum": 0.0, "maximum": 1.0},
                "length_mean": {"type": ["number", "null"], "minimum": 1},
                "length_sd": {"type": ["number", "null"], "minimum": 0},
                "length_min": {"type": ["number", "null"], "minimum": 1},
                "length_max": {"type": ["number", "null"], "minimum": 1},
            },
            "required": [
                "model_type",
                "model_file",
                "coverage",
                "pass_num",
                "min_passes",
                "min_rq",
            ],
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
                    "type": {"type": "string"},
                    "citation": {"type": "string"},
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
                "additionalProperties": False,  # Strict validation - only allow defined fields
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
                "pbsim3": {"type": "string"},
                "ccs": {"type": "string"},
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
                "coverage": {"type": "number"},
                "downsample_seed": {"type": "number"},
                "downsample_mode": {"type": "string"},
                "sample_target_bed": {"type": "string"},
                "reference_assembly": {"type": "string"},
                "vntr_region_hg19": {"type": "string"},
                "vntr_region_hg38": {"type": "string"},
                "aligner": {"type": "string", "enum": ["bwa", "minimap2"]},
                "seed": {"type": ["number", "null"]},
                "keep_intermediate_files": {"type": "boolean"},
                "vntr_capture_efficiency": {
                    "type": "object",
                    "properties": {
                        "enabled": {"type": "boolean", "default": True},
                        "penalty_factor": {
                            "type": "number",
                            "minimum": 0.1,
                            "maximum": 1.0,
                            "default": 0.375,
                        },
                        "seed": {"type": "integer", "minimum": 0, "default": 42},
                        "vntr_region": {
                            "type": "object",
                            "properties": {
                                "chr": {"type": "string"},
                                "start": {"type": "integer", "minimum": 0},
                                "end": {"type": "integer", "minimum": 1},
                                "name": {"type": "string"},
                            },
                            "required": ["chr", "start", "end", "name"],
                        },
                        "capture_bed": {"type": "string"},
                        "flanking_size": {
                            "type": "integer",
                            "minimum": 1000,
                            "maximum": 50000,
                            "default": 10000,
                        },
                        "validation": {
                            "type": "object",
                            "properties": {
                                "check_duplicates": {"type": "boolean", "default": False},
                                "report_statistics": {"type": "boolean", "default": True},
                            },
                        },
                    },
                    "additionalProperties": False,
                },
            },
            "required": [
                "human_reference",
                "threads",
            ],
            "additionalProperties": False,
        },
        "snapshot_validation": {
            "type": "object",
            "additionalProperties": {
                "type": "object",
                "properties": {
                    "description": {"type": "string"},
                    "pcr": {
                        "type": "object",
                        "properties": {
                            "forward_primer": {"type": "string"},
                            "reverse_primer": {"type": "string"},
                            "reverse_needs_rc": {"type": "boolean"},
                            "max_products": {"type": "number"},
                            "size_range": {
                                "type": "object",
                                "properties": {
                                    "min": {"type": "number"},
                                    "max": {"type": "number"},
                                },
                                "required": ["min", "max"],
                                "additionalProperties": False,
                            },
                        },
                        "required": [
                            "forward_primer",
                            "reverse_primer",
                            "reverse_needs_rc",
                            "max_products",
                            "size_range",
                        ],
                        "additionalProperties": False,
                    },
                    "digest": {
                        "type": "object",
                        "properties": {
                            "enzyme": {"type": "string"},
                            "recognition_site": {"type": "string"},
                            "expected_survivors": {"type": "string"},
                        },
                        "required": ["enzyme", "recognition_site"],
                        "additionalProperties": False,
                    },
                    "snapshot": {
                        "type": "object",
                        "properties": {
                            "primers": {
                                "type": "object",
                                "additionalProperties": {"type": "string"},
                            },
                            "fluorophore_map": {
                                "type": "object",
                                "additionalProperties": {
                                    "type": "object",
                                    "properties": {
                                        "color": {"type": "string"},
                                        "dye": {"type": "string"},
                                    },
                                    "required": ["color", "dye"],
                                    "additionalProperties": False,
                                },
                            },
                        },
                        "required": ["primers", "fluorophore_map"],
                        "additionalProperties": False,
                    },
                    "validation": {
                        "type": "object",
                        "properties": {
                            "mutant_pattern": {"type": "string"},
                            "normal_pattern": {"type": "string"},
                            "expected_mutant_fluorescence": {"type": "string"},
                        },
                        "required": [
                            "mutant_pattern",
                            "normal_pattern",
                        ],
                        "additionalProperties": False,
                    },
                },
                "required": ["pcr", "digest", "snapshot", "validation"],
                "additionalProperties": False,
            },
        },
        "toxic_protein_detection": {
            "type": "object",
            "properties": {
                "consensus_motif": {"type": "string"},
                "identity_threshold": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                "expected_repeat_count": {"type": "number", "minimum": 1},
                "key_residues": {"type": "array", "items": {"type": "string"}},
                "weights": {
                    "type": "object",
                    "properties": {
                        "repeat": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                        "composition": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                    },
                    "required": ["repeat", "composition"],
                    "additionalProperties": False,
                },
                "toxic_cutoff": {"type": "number", "minimum": 0.0, "maximum": 1.0},
            },
            "required": [
                "consensus_motif",
                "identity_threshold",
                "expected_repeat_count",
                "key_residues",
                "weights",
                "toxic_cutoff",
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
    """Load and validate MucOneUp configuration from JSON file.

    Validates configuration against CONFIG_SCHEMA and normalizes constants
    format (converts flat format to nested hg19/hg38 structure if needed).
    Performs additional validation for mutation allowed_repeats.

    Args:
        config_path: Path to the JSON config file

    Returns:
        Validated configuration dictionary containing all required sections
        (repeats, constants, probabilities, length_model, mutations,
        tools, read_simulation, nanosim_params)

    Raises:
        FileNotFoundError: If the config file does not exist
        json.JSONDecodeError: If the config file is not valid JSON
        ValidationError: If the config does not conform to CONFIG_SCHEMA or
                        contains invalid mutation allowed_repeats

    Note:
        Flat constants format (for backward compatibility) is automatically
        converted to nested format with default reference assembly (hg38).
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

        # AUTO-MIGRATION (Issue #28): Create reference_genomes from constants if missing
        if "reference_genomes" not in config and "constants" in config:
            logging.info("Auto-migrating config: adding reference_genomes section")

            config["reference_genomes"] = {}

            # Migrate from nested constants format
            for assembly in ["hg38", "hg19"]:
                if assembly in config["constants"]:
                    assembly_const = config["constants"][assembly]
                    if isinstance(assembly_const, dict) and "vntr_region" in assembly_const:
                        # Create minimal reference_genomes entry
                        config["reference_genomes"][assembly] = {
                            "fasta_path": f"reference/{assembly}/{assembly}.fa",
                            "vntr_region": assembly_const["vntr_region"],
                            "display_name": f"Auto-migrated {assembly.upper()}",
                        }
                        logging.debug(f"Auto-migrated {assembly} to reference_genomes")

            if config["reference_genomes"]:
                logging.warning(
                    "Config format is outdated. Please add 'reference_genomes' section "
                    "with proper fasta_path values. See documentation for details."
                )

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

    return config  # type: ignore[no-any-return]

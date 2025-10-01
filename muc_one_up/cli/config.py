"""
Configuration and setup functions for MucOneUp CLI.

Single Responsibility: Handle configuration loading, simulation mode determination,
and mutation configuration processing.
"""

import itertools
import json
import logging
import random
import sys
from pathlib import Path
from typing import Any

from ..io import parse_vntr_structure_file


def parse_fixed_lengths(fixed_lengths_args, num_haplotypes):
    """
    Parse the fixed-length values provided as strings.
    Each value may be a single integer or a range in the format "start-end".
    Returns a list of lists, one per haplotype.
    """
    if len(fixed_lengths_args) not in (1, num_haplotypes):
        logging.error("--fixed-lengths must have either 1 value or N=--num-haplotypes")
        sys.exit(1)
    result = []
    for val in fixed_lengths_args:
        if "-" in val:
            try:
                start_str, end_str = val.split("-")
                start = int(start_str)
                end = int(end_str)
                result.append(list(range(start, end + 1)))
            except ValueError:
                logging.error("Invalid fixed-length range format: '%s'", val)
                sys.exit(1)
        else:
            try:
                fixed = int(val)
                result.append([fixed])
            except ValueError:
                logging.error("Invalid fixed-length value: '%s'", val)
                sys.exit(1)
    if len(result) == 1 and num_haplotypes > 1:
        result = result * num_haplotypes
    return result


def build_cartesian_fixed_length_configs(fixed_matrix):
    """
    Given a list of lists (each corresponding to one haplotype's possible fixed lengths),
    return the Cartesian product as a list of fixed-length configurations (each a list).
    """
    return [list(prod) for prod in itertools.product(*fixed_matrix)]


def numbered_filename(
    out_dir: str, out_base: str, iteration: int, file_type: str, variant: str = ""
) -> str:
    """
    Build a filename by combining the output directory, base name, iteration number,
    variant suffix, and a file type suffix.
    """
    iter_str = f".{iteration:03d}"
    variant_str = f".{variant}" if variant else ""
    return str(Path(out_dir) / f"{out_base}{iter_str}{variant_str}.{file_type}")


def setup_configuration(args) -> tuple[dict[str, Any], str, str]:
    """
    Load configuration and setup output directory.

    Single Responsibility: Configuration loading and validation.

    Args:
        args: Parsed command-line arguments

    Returns:
        Tuple of (config_dict, out_dir, out_base)

    Raises:
        SystemExit: If config cannot be loaded
    """
    try:
        with Path(args.config).open() as fh:
            config = json.load(fh)
        logging.info("Configuration loaded from %s", args.config)

        if args.reference_assembly:
            current_assembly = config.get("reference_assembly", "hg38")
            config["reference_assembly"] = args.reference_assembly
            logging.info(
                "Reference assembly overridden by command line: %s -> %s",
                current_assembly,
                args.reference_assembly,
            )
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logging.error("Could not load config: %s", e)
        sys.exit(1)

    out_dir = args.out_dir
    out_base = args.out_base
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    return config, out_dir, out_base


def determine_simulation_mode(args, config) -> tuple[list, list | None, dict | None]:
    """
    Determine simulation mode and return appropriate configurations.

    Single Responsibility: Resolve simulation mode (structure file, fixed lengths, or random).

    Args:
        args: Parsed command-line arguments
        config: Configuration dictionary

    Returns:
        Tuple of (simulation_configs, predefined_chains, structure_mutation_info)
    """
    predefined_chains = None
    structure_mutation_info = None

    # Process --input-structure if provided
    if args.input_structure is not None:
        try:
            logging.info("Using input structure file: %s", args.input_structure)
            predefined_chains, structure_mutation_info = parse_vntr_structure_file(
                args.input_structure, config
            )
            num_chains = len(predefined_chains)
            logging.info("Loaded %d haplotype chains from structure file", num_chains)

            if structure_mutation_info:
                logging.info(
                    "Found mutation information in structure file: %s",
                    structure_mutation_info["name"],
                )

            if num_chains != args.num_haplotypes:
                logging.warning(
                    "--num-haplotypes=%d specified but structure file has %d chains. "
                    "Using the number of chains from the structure file.",
                    args.num_haplotypes,
                    num_chains,
                )

            if args.simulate_series is not None:
                logging.warning("--simulate-series is ignored when using --input-structure")

            simulation_configs = ["from_structure"]
            logging.info("Using predefined VNTR chains from structure file.")
        except Exception as e:
            logging.error("Error parsing structure file: %s", e)
            sys.exit(1)

    # Process --fixed-lengths if --input-structure is not provided
    elif args.fixed_lengths is not None:
        fixed_matrix = parse_fixed_lengths(args.fixed_lengths, args.num_haplotypes)
        if args.simulate_series is not None:
            step = args.simulate_series
            new_fixed_matrix = []
            for lst in fixed_matrix:
                if len(lst) > 1:
                    start = lst[0]
                    end = lst[-1]
                    new_lst = list(range(start, end + 1, step))
                    if new_lst[-1] != end:
                        new_lst.append(end)
                    new_fixed_matrix.append(new_lst)
                else:
                    new_fixed_matrix.append(lst)
            fixed_matrix = new_fixed_matrix
            simulation_configs = build_cartesian_fixed_length_configs(fixed_matrix)
            logging.info(
                "Series mode enabled with step %d: %d simulation iterations generated from fixed-length ranges.",
                step,
                len(simulation_configs),
            )
        else:
            simulation_configs = [[random.choice(lst) for lst in fixed_matrix]]
            logging.info(
                "Single simulation iteration generated using a random choice from each fixed-length range."
            )
    else:
        simulation_configs = [None]  # Use random lengths if not provided

    return simulation_configs, predefined_chains, structure_mutation_info


def process_mutation_config(
    args, structure_mutation_info
) -> tuple[bool, list[str] | None, str | None]:
    """
    Process mutation configuration from CLI args and structure file.

    Single Responsibility: Parse and validate mutation settings.

    Args:
        args: Parsed command-line arguments
        structure_mutation_info: Mutation info from structure file (if any)

    Returns:
        Tuple of (dual_mutation_mode, mutation_pair, mutation_name)
    """
    dual_mutation_mode = False
    mutation_pair = None

    # Convert structure file mutation info to CLI-style mutation info if present
    if structure_mutation_info:
        args.mutation_name = structure_mutation_info["name"]
        target_tuples = structure_mutation_info["targets"]
        args.mutation_targets = [f"{hap_idx},{rep_idx}" for hap_idx, rep_idx in target_tuples]
        logging.info(
            "Using mutation information from structure file: %s at targets %s",
            args.mutation_name,
            args.mutation_targets,
        )

    # Process standard mutation args
    if args.mutation_name:
        if "," in args.mutation_name:
            mutation_pair = [s.strip() for s in args.mutation_name.split(",")]
            if mutation_pair[0].lower() != "normal":
                logging.error("In dual simulation mode, the first mutation-name must be 'normal'.")
                sys.exit(1)
            dual_mutation_mode = True
            logging.info("Using mutations from command line: %s", mutation_pair)
        else:
            logging.info("Using mutation from command line: %s", args.mutation_name)

    return dual_mutation_mode, mutation_pair, args.mutation_name

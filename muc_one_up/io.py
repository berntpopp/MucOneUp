#!/usr/bin/env python3
"""
io.py

This module implements input/output operations for MucOneUp.
It supports:
  - Parsing VNTR structure files for simulation input
"""

import logging
import re
from typing import List, Dict, Tuple, Optional, Any, Union
import ast


def parse_vntr_structure_file(
    filepath: str, config: dict
) -> Tuple[List[List[str]], Optional[Dict[str, Any]]]:
    """
    Parse a VNTR structure file for use in simulation.

    The file format should be:
    haplotype_1<TAB>1-2-3-4-5-C-X-X-A-...-6p-7-8-9
    haplotype_2<TAB>1-2-3-4-5-C-X-X-X-...-6-7-8-9

    Each line represents one haplotype, and the symbols should be separated by
    dashes (-). The first column is the haplotype identifier, and the second
    column is the chain of repeat symbols.

    :param filepath: Path to the VNTR structure file.
    :param config: Configuration dict with valid repeat symbols.
    :return: List of lists, where each inner list contains symbols for one haplotype.
    :raises ValueError: If a symbol in the structure file is not found in the config.
    :raises FileNotFoundError: If the structure file cannot be found or read.
    """
    valid_repeat_symbols = set(config["repeats"].keys())
    haplotype_chains = []
    mutation_info = None

    try:
        # First pass: collect all comments for mutation extraction
        mutation_comments = []
        data_lines = []

        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    mutation_comments.append(line[1:].strip())
                else:
                    data_lines.append(line)

        # Extract mutation information from comments
        mutation_info = extract_mutation_info_from_comments(mutation_comments)
        if mutation_info:
            logging.info(
                f"Found mutation information in structure file: {mutation_info['name']}"
            )

        if not data_lines:
            raise ValueError(f"No valid data found in structure file: {filepath}")

        logging.info(f"Found {len(data_lines)} haplotype entries in structure file")

        for i, line in enumerate(data_lines, start=1):
            try:
                # Split by tab and get the second column (the chain)
                columns = line.split("\t")
                if len(columns) < 2:
                    raise ValueError(
                        f"Invalid format in line {i}: Expected tab-separated columns"
                    )

                chain_str = columns[1]
                chain = chain_str.split("-")

                # Validate each symbol in the chain
                for symbol in chain:
                    # Strip any mutation marker ('m') for validation
                    # This allows structures with pre-marked mutations
                    pure_symbol = symbol.rstrip("m")
                    if pure_symbol not in valid_repeat_symbols:
                        raise ValueError(
                            f"Invalid repeat symbol '{pure_symbol}' in line {i}. "
                            f"Valid symbols are: {sorted(valid_repeat_symbols)}"
                        )

                haplotype_chains.append(chain)
                logging.debug(f"Parsed haplotype {i} with {len(chain)} repeat units")

            except Exception as e:
                logging.error(f"Error parsing line {i} in structure file: {e}")
                raise ValueError(f"Failed to parse structure file at line {i}: {e}")

        return haplotype_chains, mutation_info

    except FileNotFoundError:
        logging.error(f"Structure file not found: {filepath}")
        raise FileNotFoundError(f"Structure file not found: {filepath}")
    except Exception as e:
        logging.error(f"Error reading structure file: {e}")
        raise ValueError(f"Failed to parse structure file: {e}")


def extract_mutation_info_from_comments(
    comments: List[str],
) -> Optional[Dict[str, Any]]:
    """
    Extract mutation information from structure file comments.

    Looks for comments in the format:
    # Mutation Applied: dupC (Targets: [(1, 25)])

    :param comments: List of comment strings (without the leading '#')
    :return: Dictionary with mutation information or None if no valid info found
    """
    if not comments:
        return None

    mutation_info = None

    # Look for mutation information in the comment lines
    mutation_line_pattern = re.compile(
        r"Mutation\s+Applied:\s+([\w\d]+)\s*\(\s*Targets:\s*(.+)\)"
    )

    for comment in comments:
        match = mutation_line_pattern.search(comment)
        if match:
            mutation_name = match.group(1).strip()
            targets_str = match.group(2).strip()

            try:
                # Parse the targets string as a Python literal (list of tuples)
                targets = ast.literal_eval(targets_str)

                # Ensure the targets are in the expected format
                if not all(isinstance(t, tuple) and len(t) == 2 for t in targets):
                    logging.warning(
                        f"Malformed mutation targets in structure file: {targets_str}"
                    )
                    continue

                mutation_info = {"name": mutation_name, "targets": targets}
                logging.info(
                    f"Extracted mutation information: {mutation_name} at {targets}"
                )
                break  # Use the first valid mutation info we find

            except (ValueError, SyntaxError) as e:
                logging.warning(
                    f"Failed to parse mutation targets from structure file: {e}"
                )
                continue

    return mutation_info

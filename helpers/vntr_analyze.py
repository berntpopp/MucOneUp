#!/usr/bin/env python3
"""
Helper script for analyzing VNTR repeat units using a configuration file.

This script processes a CSV/TSV file containing VNTR structure sequences,
removing duplicate sequences, and computes statistics on the repeat units
(e.g. min, max, mean, and median lengths). It also calculates a transition
probability matrix for repeat tokens (with a transition from the last token to
an "END" state).

The script expects a configuration JSON file (provided via the --config option)
that contains at least the "repeats" key with a dictionary of known repeat units.
Other keys such as "constants", "length_model", "mutations", etc., may also be
present.

Usage:
    ./vntr_analyze.py input.tsv --config config.json --header --structure-column vntr [--output output.json]

If the input file has no header, specify the column index (0-based) via
--structure-column.
"""

import argparse
import csv
import json
import re
import statistics
import sys
import warnings
from collections import defaultdict
from typing import Any, TextIO


def parse_vntr(vntr_str: str) -> list[str]:
    """
    Parse a VNTR string into a list of repeat unit tokens.

    Splits the input string using a regular expression that matches common
    dash characters (including Unicode variants).

    Args:
        vntr_str: A string representing VNTR repeat units separated by dashes.

    Returns:
        A list of repeat unit tokens.
    """
    # Match various Unicode hyphen/dash characters (U+002D, U+2010, U+2013, U+2014)
    pattern = r"[-‐–—]"  # noqa: RUF001 - Intentionally matching multiple hyphen types
    tokens = re.split(pattern, vntr_str.strip())
    return [token for token in tokens if token]


def load_config(config_file: TextIO) -> dict[str, Any]:
    """
    Load a configuration JSON file.

    Args:
        config_file: A file-like object containing the configuration JSON.

    Returns:
        A dictionary representing the configuration.
    """
    return json.load(config_file)


def analyze_vntr_sequences(
    file_handle: TextIO,
    structure_column: str | int,
    has_header: bool,
    known_repeats: dict[str, str],
    delimiter: str = "\t",
) -> dict[str, Any]:
    """
    Analyze VNTR sequences from an input file.

    This function reads VNTR structure sequences from a CSV/TSV file, removes
    duplicates, calculates the number of repeat units per sequence, and builds a
    transition probability matrix from the tokenized repeat units.

    It also checks if each token is in the known repeats dictionary and warns
    if any unknown tokens are encountered.

    Args:
        file_handle: A file-like object for the input data.
        structure_column: If has_header is True, then a column name; otherwise
            a zero-based column index (as an integer or a string convertible to int).
        has_header: True if the input file contains a header row.
        known_repeats: Dictionary of known repeat units from the configuration.
        delimiter: Field delimiter (default is tab).

    Returns:
        A dictionary with the following keys:
          - "min_repeats": Minimum number of repeat units found.
          - "max_repeats": Maximum number of repeat units found.
          - "mean_repeats": Mean number of repeat units.
          - "median_repeats": Median number of repeat units.
          - "probabilities": A transition probability matrix (including an "END" state).
    """
    seen = set()
    lengths: list[int] = []
    transition_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    unknown_tokens = set()

    if has_header:
        reader = csv.DictReader(file_handle, delimiter=delimiter)
        for row in reader:
            structure = row.get(str(structure_column), "").strip()
            if not structure:
                continue
            if structure in seen:
                continue
            seen.add(structure)
            tokens = parse_vntr(structure)
            if not tokens:
                continue
            # Check each token against known repeats.
            for token in tokens:
                if token not in known_repeats:
                    unknown_tokens.add(token)
            lengths.append(len(tokens))
            for i in range(len(tokens) - 1):
                src = tokens[i]
                dst = tokens[i + 1]
                transition_counts[src][dst] += 1
            # Transition from the last token to the "END" state.
            transition_counts[tokens[-1]]["END"] += 1
    else:
        reader = csv.reader(file_handle, delimiter=delimiter)
        try:
            col_index = int(structure_column)
        except ValueError as e:
            raise ValueError(
                "When no header is present, --structure-column must be an integer index."
            ) from e
        for row in reader:
            if not row:
                continue
            if col_index >= len(row):
                continue
            structure = row[col_index].strip()
            if not structure:
                continue
            if structure in seen:
                continue
            seen.add(structure)
            tokens = parse_vntr(structure)
            if not tokens:
                continue
            for token in tokens:
                if token not in known_repeats:
                    unknown_tokens.add(token)
            lengths.append(len(tokens))
            for i in range(len(tokens) - 1):
                src = tokens[i]
                dst = tokens[i + 1]
                transition_counts[src][dst] += 1
            transition_counts[tokens[-1]]["END"] += 1

    if unknown_tokens:
        warnings.warn(
            "The following repeat tokens were not found in the known repeats: "
            f"{', '.join(sorted(unknown_tokens))}",
            stacklevel=2,
        )

    if not lengths:
        raise ValueError("No valid VNTR sequences found in the input file.")

    min_repeats = min(lengths)
    max_repeats = max(lengths)
    mean_repeats = statistics.mean(lengths)
    median_repeats = statistics.median(lengths)

    probabilities: dict[str, dict[str, float]] = {}
    for src, dests in transition_counts.items():
        total = sum(dests.values())
        probabilities[src] = {dst: count / total for dst, count in dests.items()}

    return {
        "min_repeats": min_repeats,
        "max_repeats": max_repeats,
        "mean_repeats": mean_repeats,
        "median_repeats": median_repeats,
        "probabilities": probabilities,
    }


def main() -> None:
    """
    Main function to parse command-line arguments, load configuration, and run analysis.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Analyze VNTR sequences from a CSV/TSV file, compute statistics and "
            "transition probabilities, and combine the results with a configuration "
            "of known repeats."
        )
    )
    parser.add_argument(
        "input_file",
        type=argparse.FileType("r"),
        help="Path to the input CSV/TSV file containing VNTR structures.",
    )
    parser.add_argument(
        "--config",
        type=argparse.FileType("r"),
        required=True,
        help="Path to the configuration JSON file.",
    )
    parser.add_argument(
        "--structure-column",
        default="vntr",
        help=(
            "Column name (if header exists) or column index (0-based) "
            "that contains the VNTR structure. Default is 'vntr'."
        ),
    )
    parser.add_argument(
        "--delimiter",
        default="\t",
        help="Field delimiter for the input file (default is tab).",
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="Specify if the input file has a header row.",
    )
    parser.add_argument(
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help=("Path to the output file. If not provided, results are printed to stdout."),
    )

    args = parser.parse_args()

    try:
        config = load_config(args.config)
    except Exception as err:
        parser.error(f"Error loading config file: {err}")

    # Retrieve known repeats from config.
    known_repeats = config.get("repeats", {})
    if not known_repeats:
        parser.error(
            "Configuration file does not contain any known repeats under the 'repeats' key."
        )

    try:
        analysis_results = analyze_vntr_sequences(
            args.input_file,
            args.structure_column,
            args.header,
            known_repeats,
            args.delimiter,
        )
    except Exception as err:
        parser.error(f"Error analyzing VNTR sequences: {err}")

    # Build output using the known repeats from the config and computed results.
    output = {
        "min_repeats": analysis_results["min_repeats"],
        "max_repeats": analysis_results["max_repeats"],
        "mean_repeats": analysis_results["mean_repeats"],
        "median_repeats": analysis_results["median_repeats"],
        "repeats": known_repeats,
        "probabilities": analysis_results["probabilities"],
    }

    # Write the result to the specified output file (or stdout).
    args.output.write(json.dumps(output, indent=2) + "\n")
    args.output.close()


if __name__ == "__main__":
    main()

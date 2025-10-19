"""VNTR structure analysis and transition probability calculation.

This module provides functions for analyzing VNTR (Variable Number Tandem Repeat)
structures from CSV/TSV files. It computes statistics on repeat unit distributions
and builds transition probability matrices showing the likelihood of one repeat
unit following another.

Functions:
    parse_vntr: Parse VNTR string into repeat unit tokens
    analyze_vntr_sequences: Analyze VNTR sequences and compute statistics

Example:
    >>> from muc_one_up.analysis.vntr_statistics import analyze_vntr_sequences
    >>> known_repeats = {"1": "ACGT...", "2": "GCTA...", "X": "ATCG..."}
    >>> with open("vntr_data.tsv") as f:
    ...     results = analyze_vntr_sequences(
    ...         f, "vntr", has_header=True, known_repeats=known_repeats
    ...     )
    >>> print(f"Mean repeats: {results['mean_repeats']}")
"""

import csv
import re
import statistics
import warnings
from collections import defaultdict
from typing import Any, TextIO


def parse_vntr(vntr_str: str) -> list[str]:
    """Parse a VNTR string into a list of repeat unit tokens.

    Splits the input string using a regular expression that matches common
    dash characters (including Unicode variants: U+002D, U+2010, U+2013, U+2014).

    Args:
        vntr_str: A string representing VNTR repeat units separated by dashes.
                  Example: "1-2-3-X-A-B-6p-7-8-9"

    Returns:
        A list of repeat unit tokens.
        Example: ["1", "2", "3", "X", "A", "B", "6p", "7", "8", "9"]

    Examples:
        >>> parse_vntr("1-2-3-X-A-B-6p-7-8-9")
        ['1', '2', '3', 'X', 'A', 'B', '6p', '7', '8', '9']

        >>> parse_vntr("  1-2-3  ")  # Handles whitespace
        ['1', '2', '3']

        >>> parse_vntr("X")  # Single token
        ['X']
    """
    # Match various Unicode hyphen/dash characters (U+002D, U+2010, U+2013, U+2014)
    pattern = r"[-‐–—]"  # noqa: RUF001 - Intentionally matching multiple hyphen types
    tokens = re.split(pattern, vntr_str.strip())
    return [token for token in tokens if token]


def analyze_vntr_sequences(
    file_handle: TextIO,
    structure_column: str | int,
    has_header: bool,
    known_repeats: dict[str, str],
    delimiter: str = "\t",
) -> dict[str, Any]:
    """Analyze VNTR sequences from an input file.

    This function reads VNTR structure sequences from a CSV/TSV file, removes
    duplicates, calculates the number of repeat units per sequence, and builds a
    transition probability matrix from the tokenized repeat units.

    The transition matrix includes an "END" state representing sequence termination.
    Unknown repeat tokens (not in known_repeats) trigger warnings but don't fail.

    Args:
        file_handle: A file-like object for the input data.
        structure_column: If has_header is True, then a column name; otherwise
            a zero-based column index (as an integer or a string convertible to int).
        has_header: True if the input file contains a header row.
        known_repeats: Dictionary of known repeat units from the configuration.
                       Example: {"1": "ACGT...", "X": "GCTA...", ...}
        delimiter: Field delimiter (default is tab).

    Returns:
        A dictionary with the following keys:
          - "min_repeats": Minimum number of repeat units found (int)
          - "max_repeats": Maximum number of repeat units found (int)
          - "mean_repeats": Mean number of repeat units (float)
          - "median_repeats": Median number of repeat units (float)
          - "probabilities": Transition probability matrix (dict[str, dict[str, float]])
                             Example: {"X": {"A": 0.5, "B": 0.3, "END": 0.2}}

    Raises:
        ValueError: If no valid VNTR sequences found in input file.
        ValueError: If structure_column is invalid when has_header=False.

    Warns:
        UserWarning: If unknown repeat tokens encountered (not in known_repeats).

    Examples:
        >>> with open("vntr_data.tsv") as f:
        ...     results = analyze_vntr_sequences(
        ...         f, "vntr", True, {"1": "ACGT", "2": "GCTA"}, "\t"
        ...     )
        >>> print(results["mean_repeats"])
        42.5

        >>> # Build transition matrix
        >>> probs = results["probabilities"]
        >>> print(probs["X"]["A"])  # Probability of X -> A transition
        0.65
    """
    seen = set()
    lengths: list[int] = []
    transition_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    unknown_tokens = set()

    if has_header:
        dict_reader = csv.DictReader(file_handle, delimiter=delimiter)
        for row in dict_reader:
            structure = row.get(str(structure_column), "").strip()
            if not structure:
                continue
            if structure in seen:
                continue
            seen.add(structure)
            tokens = parse_vntr(structure)
            if not tokens:
                continue
            # Check each token against known repeats
            for token in tokens:
                if token not in known_repeats:
                    unknown_tokens.add(token)
            lengths.append(len(tokens))
            # Build transition matrix
            for i in range(len(tokens) - 1):
                src = tokens[i]
                dst = tokens[i + 1]
                transition_counts[src][dst] += 1
            # Transition from the last token to the "END" state
            transition_counts[tokens[-1]]["END"] += 1
    else:
        csv_reader = csv.reader(file_handle, delimiter=delimiter)
        try:
            col_index = int(structure_column)
        except ValueError as e:
            raise ValueError(
                "When no header is present, --structure-column must be an integer index."
            ) from e
        for csv_row in csv_reader:
            if not csv_row:
                continue
            if col_index >= len(csv_row):
                continue
            structure = csv_row[col_index].strip()
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

    # Warn about unknown tokens
    if unknown_tokens:
        warnings.warn(
            "The following repeat tokens were not found in the known repeats: "
            f"{', '.join(sorted(unknown_tokens))}",
            stacklevel=2,
        )

    if not lengths:
        raise ValueError("No valid VNTR sequences found in the input file.")

    # Calculate statistics
    min_repeats = min(lengths)
    max_repeats = max(lengths)
    mean_repeats = statistics.mean(lengths)
    median_repeats = statistics.median(lengths)

    # Convert counts to probabilities
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

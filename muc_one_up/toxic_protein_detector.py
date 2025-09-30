#!/usr/bin/env python3
"""
toxic_protein_detector.py

This module implements detection of a toxic protein sequence from ORF FASTA outputs.
The detection algorithm is based on quantifying how far the protein`s variable (repeat)
region deviates from the wild-type. The analysis proceeds through four main steps:

1. **Preprocessing & Region Extraction:**
   - If constant left/right flanks are provided, they are removed from the ORF sequence,
     leaving only the variable (repeat) region for analysis.
   - If the mutated region is known (or determined by alignment with a wild-type),
     it is extracted accordingly.

2. **Detection and Quantification of Repeats:**
   - A sliding window (with window length equal to the consensus repeat motif) is applied
     to the variable region.
   - For each window, the similarity (using a simple Hamming distance measure) to the
     consensus motif (e.g. "RCHLGPGHQAGPGLHR") is computed.
   - Windows with similarity above a defined threshold (default 80% identity) are considered
     valid repeat matches.
   - The number of such matches (repeat_count) and their average identity (avg_repeat_identity)
     are recorded.
   - A repeat score is calculated by scaling the average identity by the expected repeat count.

3. **Amino Acid Composition Analysis:**
   - The frequencies of key amino acids (by default R, C, and H) in the variable region are computed.
   - A wild-type composition is defined (by default, the consensus motif repeated a given number of times)
     and the similarity score is computed as:

       S_composition = 1 - (sum(|f_mut - f_wt|) / sum(f_wt))

     where a score near 1 indicates high similarity in composition.

4. **Combining Metrics into an Overall Score:**
   - The overall deviation (or similarity) score is calculated as a weighted sum of the repeat score and
     the composition similarity. In this implementation, a **higher overall score indicates a greater deviation
     from the wild-type (i.e. the protein is toxic)**.
   - If the overall score exceeds a user-defined toxic cutoff (e.g. 0.5), the ORF is flagged as toxic.

The module also provides a function to scan an ORF FASTA file and generate a dictionary mapping each ORF header
to its detection metrics.

Usage (via command line):
    python toxic_protein_detector.py <orf_fasta> [options]

The output is a JSON-formatted report containing, for each ORF:
    - repeat_count
    - avg_repeat_identity
    - repeat_score
    - composition_similarity
    - overall_score
    - toxic_flag (1.0 if flagged as toxic, 0.0 otherwise)

All parameters (consensus motif, identity threshold, expected repeat count, weights, and cutoff) are configurable.
"""

import argparse
import json
from pathlib import Path


def sliding_window_repeat_analysis(
    region: str, consensus: str, identity_threshold: float = 0.8
) -> tuple[int, float]:
    """
    Slide a window of length equal to the consensus motif along the given region
    and compute repeat matches based on identity.

    For each window (substring) of the protein variable region, the function computes the
    fraction of positions that match the consensus motif. If the identity is above the
    specified threshold, the window is recorded as a valid repeat match.

    Args:
        region (str): Protein sequence (variable region) to analyze.
        consensus (str): Consensus motif sequence (e.g., "RCHLGPGHQAGPGLHR").
        identity_threshold (float): Minimum identity (0 to 1) to count a window as a match.
                                    Default is 0.8 (80% identity).

    Returns:
        Tuple[int, float]: A tuple containing:
            - The number of detected repeat matches.
            - The average identity of all detected repeats.
    """
    motif_len = len(consensus)
    detected_repeats = []
    for i in range(len(region) - motif_len + 1):
        candidate = region[i : i + motif_len]
        identity = sum(1 for a, b in zip(candidate, consensus, strict=True) if a == b) / motif_len
        if identity >= identity_threshold:
            detected_repeats.append((i, candidate, identity))
    repeat_count = len(detected_repeats)
    avg_identity = (
        sum(rep[2] for rep in detected_repeats) / repeat_count if repeat_count > 0 else 0.0
    )
    return repeat_count, avg_identity


def compute_amino_acid_composition(seq: str, residues: list[str]) -> dict[str, float]:
    """
    Compute the frequency of each specified amino acid in a protein sequence.

    The frequency is computed as the count of the amino acid divided by the total length
    of the sequence.

    Args:
        seq (str): Protein sequence.
        residues (List[str]): List of amino acids to include in the analysis.

    Returns:
        Dict[str, float]: A dictionary mapping each amino acid to its frequency in the sequence.
    """
    total = len(seq)
    freq = {}
    for aa in residues:
        freq[aa] = seq.count(aa) / total if total > 0 else 0.0
    return freq


def composition_similarity(mutant_seq: str, wildtype_seq: str, key_residues: list[str]) -> float:
    """
    Compute a similarity score based on the frequencies of key amino acids
    between the mutant (toxic) and wild-type protein sequences.

    The score is calculated as:
        1 - (sum(|f_mut - f_wt|) / sum(f_wt))
    so that a score close to 1 indicates that the composition is very similar to wild-type,
    while lower scores indicate divergence.

    Args:
        mutant_seq (str): Protein sequence from the toxic variant.
        wildtype_seq (str): Protein sequence from the wild-type.
        key_residues (List[str]): List of key amino acids (e.g., ['R', 'C', 'H']).

    Returns:
        float: Composition similarity score (between 0 and 1).
    """
    comp_mut = compute_amino_acid_composition(mutant_seq, key_residues)
    comp_wt = compute_amino_acid_composition(wildtype_seq, key_residues)
    diff = sum(abs(comp_mut[aa] - comp_wt[aa]) for aa in key_residues)
    norm = sum(comp_wt[aa] for aa in key_residues)
    return 1 - diff / norm if norm > 0 else 0.0


def detect_toxic_protein_in_sequence(
    protein_seq: str,
    left_const: str | None = None,
    right_const: str | None = None,
    consensus: str = "RCHLGPGHQAGPGLHR",
    identity_threshold: float = 0.8,
    key_residues: list[str] | None = None,
    expected_repeat_count: int = 10,
    w_repeat: float = 0.6,
    w_composition: float = 0.4,
    wildtype_variable_region: str | None = None,
    toxic_detection_cutoff: float = 0.5,
) -> dict[str, float]:
    """
    Analyze a protein sequence for toxic protein features by quantifying its repeat structure
    and amino acid composition relative to a wild-type model.

    **Step 1: Region Extraction**
      - If constant left/right flanks are provided, they are removed from the sequence,
        leaving only the variable (repeat) region for analysis.

    **Step 2: Repeat Analysis**
      - A sliding window (of length equal to the consensus motif) is applied to the variable region.
      - The identity of each window compared to the consensus motif is computed.
      - Windows with identity above the threshold are considered valid repeats.
      - The number of repeats and the average identity are computed and combined into a repeat score,
        scaled by the expected number of repeats.

    **Step 3: Composition Analysis**
      - The frequency of key residues (default: R, C, H) is computed for both the variable region and
        a wild-type model (by default, the consensus motif repeated).
      - A composition similarity score is computed such that values near 1 indicate high similarity to
        the wild-type.

    **Step 4: Combining Metrics**
      - The overall deviation score is calculated as a weighted sum of the repeat score and the composition
        similarity score.
      - In this implementation, a higher overall score indicates a greater deviation from the wild-type,
        and thus flags the protein as toxic.
      - If the overall score exceeds the specified toxic detection cutoff, the function flags the protein
        as toxic (toxic_flag = 1.0); otherwise, it is considered normal (toxic_flag = 0.0).

    Args:
        protein_seq (str): Full protein sequence from an ORF.
        left_const (Optional[str]): Left constant region to remove (if provided).
        right_const (Optional[str]): Right constant region to remove (if provided).
        consensus (str): Consensus motif used for repeat detection.
        identity_threshold (float): Identity threshold for a window to be considered a repeat.
        key_residues (Optional[List[str]]): List of key amino acids for composition analysis.
                                          Defaults to ['R', 'C', 'H'] if not provided.
        expected_repeat_count (int): Expected number of repeats in the wild-type.
        w_repeat (float): Weight for the repeat analysis score.
        w_composition (float): Weight for the composition similarity score.
        wildtype_variable_region (Optional[str]): Wild-type variable region sequence. If not provided,
                                                  it is approximated as the consensus motif repeated.
        toxic_detection_cutoff (float): Overall score cutoff above which the protein is flagged as toxic.

    Returns:
        Dict[str, float]: A dictionary containing:
            - 'repeat_count': Detected number of repeat matches.
            - 'avg_repeat_identity': Average identity of the detected repeats.
            - 'repeat_score': Score derived from repeat analysis.
            - 'composition_similarity': Amino acid composition similarity score.
            - 'overall_score': Combined overall deviation score.
            - 'toxic_flag': 1.0 if the overall score exceeds the cutoff (toxic), 0.0 otherwise.
    """
    if key_residues is None:
        key_residues = ["R", "C", "H"]

    # Extract the variable region by removing constant flanks if provided.
    variable_region = protein_seq
    if left_const and protein_seq.startswith(left_const):
        variable_region = protein_seq[len(left_const) :]
    if right_const and variable_region.endswith(right_const):
        variable_region = variable_region[: -len(right_const)]

    # --- Step 2: Repeat Analysis ---
    # Slide a window over the variable region and compute the identity with the consensus motif.
    repeat_count, avg_identity = sliding_window_repeat_analysis(
        variable_region, consensus, identity_threshold
    )
    # Scale the average identity by the minimum of detected repeats and expected repeat count.
    repeat_score = avg_identity * min(repeat_count, expected_repeat_count) / expected_repeat_count

    # --- Step 3: Composition Analysis ---
    # Define the wild-type variable region if not provided.
    if wildtype_variable_region is None:
        wildtype_variable_region = consensus * expected_repeat_count
    # Compute the composition similarity between the variable region and the wild-type.
    comp_sim = composition_similarity(variable_region, wildtype_variable_region, key_residues)

    # --- Step 4: Combine Metrics ---
    overall_score = w_repeat * repeat_score + w_composition * comp_sim

    # Flag the protein as toxic if the overall score exceeds the toxic_detection_cutoff.
    toxic_flag = 1.0 if overall_score > toxic_detection_cutoff else 0.0

    return {
        "repeat_count": repeat_count,
        "avg_repeat_identity": avg_identity,
        "repeat_score": repeat_score,
        "composition_similarity": comp_sim,
        "overall_score": overall_score,
        "toxic_flag": toxic_flag,
    }


def scan_orf_fasta(
    orf_fasta_path: str,
    left_const: str | None = None,
    right_const: str | None = None,
    **detection_kwargs,
) -> dict[str, dict[str, float]]:
    """
    Scan an ORF FASTA file for toxic protein sequence features.

    For each ORF in the FASTA file, this function:
      - Reads the sequence (ignoring header lines starting with '>').
      - Calls `detect_toxic_protein_in_sequence()` to compute the detection metrics.
      - Returns a dictionary mapping each ORF header to its corresponding metrics.

    Args:
        orf_fasta_path (str): Path to the ORF FASTA file.
        left_const (Optional[str]): Left constant region to remove (if applicable).
        right_const (Optional[str]): Right constant region to remove (if applicable).
        **detection_kwargs: Additional keyword arguments to be passed to detect_toxic_protein_in_sequence().

    Returns:
        Dict[str, Dict[str, float]]: A dictionary where each key is an ORF header and each value is
                                      a dictionary of detection metrics.
    """
    results = {}
    header = None
    seq_lines = []
    with Path(orf_fasta_path).open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_lines:
                    protein_seq = "".join(seq_lines)
                    metrics = detect_toxic_protein_in_sequence(
                        protein_seq, left_const, right_const, **detection_kwargs
                    )
                    results[header] = metrics
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:
            protein_seq = "".join(seq_lines)
            metrics = detect_toxic_protein_in_sequence(
                protein_seq, left_const, right_const, **detection_kwargs
            )
            results[header] = metrics
    return results


def main() -> None:
    """
    Command-line interface for toxic protein detection.

    This function parses command-line arguments, calls scan_orf_fasta() with the provided options,
    and prints the detection results as a JSON-formatted string.

    Command-line options include:
      - Path to the ORF FASTA file.
      - Optional left/right constant regions.
      - Consensus motif and identity threshold for repeat detection.
      - Expected repeat count and toxic cutoff for flagging.
    """
    parser = argparse.ArgumentParser(
        description="Scan an ORF FASTA file for toxic protein sequence features."
    )
    parser.add_argument("orf_fasta", help="Path to the ORF FASTA file.")
    parser.add_argument("--left-const", help="Left constant region (to remove).", default=None)
    parser.add_argument("--right-const", help="Right constant region (to remove).", default=None)
    parser.add_argument(
        "--consensus",
        help="Consensus motif for repeat detection.",
        default="RCHLGPGHQAGPGLHR",
    )
    parser.add_argument(
        "--identity-threshold",
        type=float,
        help="Identity threshold for repeat matching.",
        default=0.8,
    )
    parser.add_argument(
        "--expected-repeat-count",
        type=int,
        help="Expected repeat count in wild-type.",
        default=10,
    )
    parser.add_argument(
        "--toxic-cutoff",
        type=float,
        help="Overall score cutoff to flag toxic protein (default=0.5).",
        default=0.5,
    )
    args = parser.parse_args()

    detection_results = scan_orf_fasta(
        args.orf_fasta,
        left_const=args.left_const,
        right_const=args.right_const,
        consensus=args.consensus,
        identity_threshold=args.identity_threshold,
        expected_repeat_count=args.expected_repeat_count,
        toxic_detection_cutoff=args.toxic_cutoff,
    )
    print(json.dumps(detection_results, indent=4))


if __name__ == "__main__":
    main()

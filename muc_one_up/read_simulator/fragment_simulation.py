#!/usr/bin/env python3
"""
Fragment simulation module.

This module contains the ported w-Wessim2 logic for simulating read fragments
from a FASTA file with systematic errors. It handles fragment picking based on
PSL alignments and produces fragment pairs for later sequencing.
"""

import logging
import random
from pathlib import Path

# Import exceptions
from ..exceptions import FileOperationError

# Map for DNA complementation
COMP_MAP = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "n": "n",
}


def read_fasta_to_dict(fasta_file: str) -> dict[str, str]:
    """
    Read a FASTA file into a dictionary mapping chromosome names to sequences.

    Args:
        fasta_file: Path to the FASTA file.

    Returns:
        Dictionary mapping chromosome names to sequences.

    Raises:
        IOError: If the file can't be read.
        ValueError: If the FASTA format is invalid.
    """
    sequences = {}
    current_chrom: str | None = None
    current_seq: list[str] = []

    try:
        with Path(fasta_file).open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # Save the previous chromosome sequence if it exists
                    if current_chrom:
                        sequences[current_chrom] = "".join(current_seq)

                    # Start new chromosome
                    current_chrom = line[1:].split()[0]  # Extract chromosome name
                    current_seq = []
                else:
                    current_seq.append(line)

            # Save the last chromosome
            if current_chrom is not None:
                sequences[current_chrom] = "".join(current_seq)
    except Exception as e:
        logging.error(f"Error reading FASTA file {fasta_file}: {e!s}")
        raise

    return sequences


def read_syser_file(
    syser_file: str,
) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    """
    Read the systematic errors file and return four dictionaries.

    Args:
        syser_file: Path to the syser file.

    Returns:
        Tuple (tendendict_f, tendendict_r, ratesdict_f, ratesdict_r).

    Raises:
        IOError: If the file can't be read.
        ValueError: If the format is invalid.
    """
    # Dictionaries to store forward/reverse tendencies and rates
    tendendict_f = {}
    tendendict_r = {}
    ratesdict_f = {}
    ratesdict_r = {}
    direction = None
    rat = False

    try:
        with Path(syser_file).open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("@") and len(line) < 500:
                    # Extract identifier from the header line
                    space_index = line.find(" ")
                    identifier = line[1:space_index] if space_index != -1 else line[1:]

                    # Determine read direction
                    if "forward" in line:
                        direction = "f"
                    elif "reverse" in line:
                        direction = "r"
                    else:
                        logging.warning("Direction not found in syser header")
                    continue
                elif line == "+":
                    rat = True
                    continue
                else:
                    if rat:
                        # Rate information
                        if direction == "f":
                            ratesdict_f[identifier] = line
                        elif direction == "r":
                            # Reverse the string for reverse reads
                            ratesdict_r[identifier] = line[::-1]
                        rat = False
                    else:
                        # Tendency information
                        if direction == "f":
                            tendendict_f[identifier] = line
                        elif direction == "r":
                            # Reverse the string for reverse reads
                            tendendict_r[identifier] = line[::-1]
    except Exception as e:
        logging.error(f"Error reading systematic errors file {syser_file}: {e!s}")
        raise

    return tendendict_f, tendendict_r, ratesdict_f, ratesdict_r


def read_psl_file(psl_file: str) -> list[tuple[str, int, int, str]]:
    """
    Read a PSL file (skipping header lines) and return a list of matches.

    Each match is a tuple: (chrom, start, end, strand).

    Args:
        psl_file: Path to the PSL file.

    Returns:
        List of match tuples.

    Raises:
        IOError: If the file can't be read.
    """
    matches = []

    try:
        with Path(psl_file).open() as f:
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()

                # Skip header lines
                if (
                    line.startswith("psLayout")
                    or line.startswith("match")
                    or line.startswith("------")
                    or not line
                ):
                    continue

                cols = line.split("\t")
                if len(cols) < 21:
                    logging.warning(f"Skipping invalid PSL line {line_num}: insufficient columns")
                    continue

                # Example PSL format: matches misMatches repMatches nCount
                # qName qSize qStart qEnd strand tName tSize tStart tEnd blocks blockSizes qStarts tStarts

                # Extract target (reference) information
                chrom = cols[13]
                start = int(cols[15])
                end = int(cols[16])
                strand = "+" if cols[8] == "+" else "-"

                matches.append((chrom, start, end, strand))
    except Exception as e:
        logging.error(f"Error reading PSL file {psl_file}: {e!s}")
        raise

    if not matches:
        logging.warning(f"No matches found in PSL file: {psl_file}")

    return matches


def pick_on_match(
    matches: list[tuple[str, int, int, str]],
) -> tuple[str, int, int, str]:
    """
    Randomly pick one match from the provided list.

    Args:
        matches: List of match tuples.

    Returns:
        A single match tuple.

    Raises:
        ValueError: If matches list is empty.
    """
    if not matches:
        raise ValueError("No matches provided to pick from")

    return random.choice(matches)


def get_insert_length(mu: float, sigma: float, lower: int) -> int:
    """
    Return an insert length sampled from a Gaussian distribution (>= lower).

    Args:
        mu: Mean insert length.
        sigma: Standard deviation.
        lower: Minimum allowed length.

    Returns:
        Integer insert length.
    """
    while True:
        length = int(random.gauss(mu, sigma))
        if length >= lower:
            return length


def pick_fragment(
    match: tuple[str, int, int, str], ins: int, bind: float
) -> tuple[str, int, int, str]:
    """
    Randomly pick a fragment based on a match and desired insert length.

    Args:
        match: Tuple (chrom, probe_start, probe_end, strand).
        ins: Desired insert length.
        bind: Minimum fraction (%) for overlap.

    Returns:
        Tuple (chrom, fragment_start, fragment_end, strand).
    """
    chrom, mstart, mend, strand = match

    # Calculate probe length and minimum base overlap
    plen = mend - mstart
    min_overlap = int(plen * bind)

    # Calculate available space for fragment
    space = ins + plen - min_overlap

    # Pick a random offset within the available space
    offset = random.randint(0, space)

    # Calculate fragment boundaries
    f_start = mstart - offset
    if f_start < 0:
        f_start = 0
    f_end = f_start + ins

    return (chrom, f_start, f_end, strand)


def comp(sequence: str) -> str:
    """
    Return the complement of a DNA sequence (preserving case).

    Args:
        sequence: DNA sequence.

    Returns:
        Complemented sequence.
    """
    comp_seq = []
    for base in sequence:
        comp_base = COMP_MAP.get(base, base)
        comp_seq.append(comp_base)

    return "".join(comp_seq)


def simulate_fragments(
    ref_fa: str,
    syser_file: str,
    psl_file: str,
    read_number: int,
    fragment_size: int,
    fragment_sd: int,
    min_fragment: int,
    bind: float,
    output_fragments: str,
) -> None:
    """
    Simulate fragments (port of w-Wessim2) and write paired fragment sequences to a FASTA file.

    Args:
        ref_fa: Reference FASTA file.
        syser_file: Systematic errors file.
        psl_file: PSL file.
        read_number: Number of read pairs to generate.
        fragment_size: Mean fragment size.
        fragment_sd: Standard deviation of fragment size.
        min_fragment: Minimum fragment size.
        bind: Minimum fraction (%) for overlap.
        output_fragments: Output FASTA filename for fragments.

    Raises:
        SystemExit: If the simulation fails or produces invalid output.
    """
    logging.info("Starting fragment simulation (ported w-Wessim2 logic)...")

    # Load reference sequences
    try:
        logging.info(f"Reading reference FASTA: {ref_fa}")
        genome = read_fasta_to_dict(ref_fa)
        logging.info(f"Loaded {len(genome)} sequences from reference")
    except Exception as e:
        logging.error(f"Failed to load reference FASTA: {e!s}")
        raise FileOperationError(f"Failed to load reference FASTA {ref_fa}: {e!s}") from e

    # Load systematic errors
    try:
        logging.info(f"Reading systematic errors file: {syser_file}")
        tendendict_f, tendendict_r, ratesdict_f, ratesdict_r = read_syser_file(syser_file)
        logging.info(
            f"Loaded {len(tendendict_f)} forward and {len(tendendict_r)} reverse positions"
        )
    except Exception as e:
        logging.error(f"Failed to load systematic errors: {e!s}")
        raise FileOperationError(
            f"Failed to load systematic errors file {syser_file}: {e!s}"
        ) from e

    # Load PSL matches
    try:
        logging.info(f"Reading PSL alignments: {psl_file}")
        matches = read_psl_file(psl_file)
        logging.info(f"Loaded {len(matches)} matches from PSL file")
    except Exception as e:
        logging.error(f"Failed to load PSL matches: {e!s}")
        raise FileOperationError(f"Failed to load PSL alignments file {psl_file}: {e!s}") from e

    if not matches:
        logging.error("No matches found in PSL file. Cannot simulate fragments.")
        raise FileOperationError(
            f"No matches found in PSL file {psl_file}. Cannot simulate fragments."
        )

    # Simulate fragments and write to FASTA file
    try:
        with Path(output_fragments).open("w") as fout:
            for i in range(read_number):
                # Pick a random match from the PSL file
                match = pick_on_match(matches)

                # Determine fragment length from Gaussian distribution
                ins = get_insert_length(fragment_size, fragment_sd, min_fragment)

                # Pick a fragment location
                fragment = pick_fragment(match, ins, bind)

                chrom, fstart, fend, strand = fragment

                # Get fragment sequence
                if chrom not in genome:
                    logging.warning(
                        f"Chromosome {chrom} not found in reference. Skipping fragment."
                    )
                    continue

                chrom_seq = genome[chrom]

                # Make sure fragment end doesn't exceed chromosome length
                if fend > len(chrom_seq):
                    fend = len(chrom_seq)

                # Skip if fragment is too short after adjustment
                if fend - fstart < min_fragment:
                    logging.debug(f"Fragment too short after adjustment: {fend - fstart} bp")
                    continue

                # Get fragment sequence
                frag_seq = chrom_seq[fstart:fend]
                frag_seq_rc = comp(frag_seq[::-1])  # Reverse complement
                frag_len = len(frag_seq)

                # Get systematic error profiles for this fragment
                # Use chromosome as key for tendency and rates dictionaries
                tendenseq_f = tendendict_f.get(chrom, "N" * frag_len)
                ratesseq_f = ratesdict_f.get(chrom, "!" * frag_len)
                tendenseq_r = tendendict_r.get(chrom, "N" * frag_len)
                ratesseq_r = ratesdict_r.get(chrom, "!" * frag_len)

                # Fix field lengths to match fragment length
                def fix_field(field, length, pad_char):
                    if len(field) > length:
                        return field[:length]
                    else:
                        return field.ljust(length, pad_char)

                tendenseq_f_fixed = fix_field(tendenseq_f, frag_len, "N")
                ratesseq_f_fixed = fix_field(ratesseq_f, frag_len, "!")
                tendenseq_r_fixed = fix_field(tendenseq_r, frag_len, "N")
                ratesseq_r_fixed = fix_field(ratesseq_r, frag_len, "!")

                # Write fragment read pair to FASTA with appropriate headers
                # reseq expects format: >id 1;length;tendencies;rates
                if strand == "+":
                    fout.write(f">{i+1} 1;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n")
                    fout.write(f"{frag_seq}\n")
                    fout.write(f">{i+1} 2;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n")
                    fout.write(f"{frag_seq_rc}\n")
                elif strand == "-":
                    fout.write(f">{i+1} 1;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n")
                    fout.write(f"{frag_seq_rc}\n")
                    fout.write(f">{i+1} 2;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n")
                    fout.write(f"{frag_seq}\n")

                # Progress logging
                if (i + 1) % 10000 == 0:
                    logging.info(f"Generated {i + 1}/{read_number} fragment pairs")

        logging.info(f"Fragment simulation completed. Output: {output_fragments}")
    except Exception as e:
        logging.error(f"Error simulating fragments: {e!s}")
        raise FileOperationError(f"Error writing fragments to {output_fragments}: {e!s}") from e

    # Check output file exists and is non-empty
    output_path = Path(output_fragments)
    if not output_path.exists() or output_path.stat().st_size == 0:
        logging.error(
            f"Fragment simulation failed: Output file missing or empty: {output_fragments}"
        )
        raise FileOperationError(
            f"Fragment simulation failed: Output file missing or empty: {output_fragments}"
        )

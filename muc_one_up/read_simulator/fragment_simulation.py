#!/usr/bin/env python3
"""
Fragment simulation module.

This module contains Wessim2-inspired fragment simulation logic, ported and
adapted from the Wessim2 approach (not calling the Wessim2 binary). It simulates
read fragments from a FASTA file with systematic errors, handling fragment
picking based on PSL alignments and producing fragment pairs for later
sequencing via ReSeq.
"""

import logging
import math
import random
from pathlib import Path

# Import exceptions
from ..exceptions import FileOperationError
from .utils.sequence_utils import reverse_complement


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
    rng: random.Random | None = None,
) -> tuple[str, int, int, str]:
    """
    Randomly pick one match from the provided list.

    Args:
        matches: List of match tuples.
        rng: Optional Random instance for reproducibility.

    Returns:
        A single match tuple.

    Raises:
        ValueError: If matches list is empty.
    """
    if not matches:
        raise ValueError("No matches provided to pick from")

    return (rng or random).choice(matches)


def get_insert_length(mu: float, sigma: float, lower: int, rng: random.Random | None = None) -> int:
    """
    Return an insert length sampled from a Gaussian distribution (>= lower).

    Args:
        mu: Mean insert length.
        sigma: Standard deviation.
        lower: Minimum allowed length.
        rng: Optional Random instance for reproducibility.

    Returns:
        Integer insert length.
    """
    if sigma <= 0 and int(mu) < lower:
        raise ValueError(
            f"Mean fragment size ({mu}) is smaller than minimum required length ({lower}) with sigma={sigma}"
        )
    _rng = rng or random
    for _ in range(10_000):
        length = int(_rng.gauss(mu, sigma))
        if length >= lower:
            return length
    raise ValueError(
        f"Could not sample insert length >= {lower} from Gaussian(mu={mu}, sigma={sigma}) after 10,000 attempts"
    )


def pick_fragment(
    match: tuple[str, int, int, str],
    ins: int,
    bind: float,
    rng: random.Random | None = None,
) -> tuple[str, int, int, str]:
    """
    Randomly pick a fragment based on a match and desired insert length.

    Args:
        match: Tuple (chrom, probe_start, probe_end, strand).
        ins: Desired insert length.
        bind: Minimum fraction (%) for overlap.
        rng: Optional Random instance for reproducibility.

    Returns:
        Tuple (chrom, fragment_start, fragment_end, strand).
    """
    if not (0.0 <= bind <= 1.0):
        raise ValueError(f"bind must be between 0.0 and 1.0, got {bind}")

    chrom, mstart, mend, strand = match

    plen = mend - mstart
    min_overlap = math.ceil(plen * bind)
    if plen <= 0 or ins < min_overlap:
        raise ValueError(
            f"Insert size ({ins}) is smaller than minimum required probe overlap ({min_overlap})"
        )

    lower = max(0, mstart + min_overlap - ins)
    upper = mend - min_overlap
    if lower > upper:
        raise ValueError(
            f"Invalid fragment bounds [{lower}, {upper}] for match {match}, ins={ins}, bind={bind}"
        )
    f_start = (rng or random).randint(lower, upper)
    f_end = f_start + ins

    return (chrom, f_start, f_end, strand)


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
    seed: int | None = None,
    fragment_origins_path: str | None = None,
) -> None:
    """
    Simulate fragments (Wessim2-inspired, ported logic) and write paired fragment sequences to a FASTA file.

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
        seed: Random seed for reproducibility. If None, uses system randomness.

    Raises:
        SystemExit: If the simulation fails or produces invalid output.
    """
    if not (0.0 <= bind <= 1.0):
        raise FileOperationError(f"Invalid bind fraction: {bind}")

    # Create local RNG for reproducibility if seed provided
    rng = random.Random(seed) if seed is not None else None
    if seed is not None:
        logging.info(f"Fragment simulation using random seed: {seed}")

    logging.info("Starting fragment simulation (Wessim2-style, ported logic)...")

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

    # Collect fragment origins for source tracking sidecar
    fragment_origins: list[dict] = []

    # Simulate fragments and write to FASTA file
    try:
        with Path(output_fragments).open("w") as fout:
            for i in range(read_number):
                # Pick a random match from the PSL file
                match = pick_on_match(matches, rng=rng)

                # Determine fragment length from Gaussian distribution
                min_overlap = math.ceil((match[2] - match[1]) * bind)
                ins = get_insert_length(
                    fragment_size, fragment_sd, max(min_fragment, min_overlap), rng=rng
                )

                # Pick a fragment location
                fragment = pick_fragment(match, ins, bind, rng=rng)

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
                frag_seq_rc = reverse_complement(frag_seq)
                frag_len = len(frag_seq)

                # Get systematic error profiles for this fragment sliced by coordinates
                # Note: read_syser_file reverses reverse-strand strings on load, so genomic
                # [fstart:fend] maps to [len - fend : len - fstart] in tendendict_r / ratesdict_r.
                tendenseq_f_fixed = tendendict_f.get(chrom, "")[fstart:fend].ljust(frag_len, "N")
                ratesseq_f_fixed = ratesdict_f.get(chrom, "")[fstart:fend].ljust(frag_len, "!")
                r_tenden = tendendict_r.get(chrom, "")
                r_rates = ratesdict_r.get(chrom, "")
                tendenseq_r_fixed = r_tenden[
                    max(0, len(r_tenden) - fend) : max(0, len(r_tenden) - fstart)
                ].ljust(frag_len, "N")
                ratesseq_r_fixed = r_rates[
                    max(0, len(r_rates) - fend) : max(0, len(r_rates) - fstart)
                ].ljust(frag_len, "!")

                # Write fragment read pair to FASTA with appropriate headers
                # reseq expects format: >id 1;length;tendencies;rates
                if strand == "+":
                    fout.write(f">{i + 1} 1;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n")
                    fout.write(f"{frag_seq}\n")
                    fout.write(f">{i + 1} 2;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n")
                    fout.write(f"{frag_seq_rc}\n")
                elif strand == "-":
                    fout.write(f">{i + 1} 1;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n")
                    fout.write(f"{frag_seq_rc}\n")
                    fout.write(f">{i + 1} 2;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n")
                    fout.write(f"{frag_seq}\n")

                # Track fragment origin for source tracking sidecar
                if fragment_origins_path is not None:
                    fragment_origins.append(
                        {
                            "fragment_index": i + 1,
                            "chrom": chrom,
                            "fstart": fstart,
                            "fend": fend,
                            "strand": strand,
                        }
                    )

                # Progress logging
                if (i + 1) % 10000 == 0:
                    logging.info(f"Generated {i + 1}/{read_number} fragment pairs")

        logging.info(f"Fragment simulation completed. Output: {output_fragments}")

        # Write fragment origins sidecar if requested
        if fragment_origins_path is not None and fragment_origins:
            from .parsers.illumina_parser import write_fragment_origins

            write_fragment_origins(fragment_origins, fragment_origins_path)
            logging.info("Fragment origins sidecar written: %s", fragment_origins_path)

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

# muc_one_up/fasta_writer.py

import logging
from pathlib import Path


def write_fasta(sequences, filename, prefix="haplotype", comment=None, comments=None):
    """
    Write a list of sequences to a FASTA file.

    :param sequences: List of sequence strings.
    :param filename: Output filename for the FASTA file.
    :param prefix: Header prefix for each sequence (default "haplotype").
    :param comment: Optional comment to append to all header lines.
    :param comments: Optional list of comments, one per sequence.
    """
    try:
        with Path(filename).open("w") as fh:
            for i, seq in enumerate(sequences, start=1):
                header = f">{prefix}_{i}"

                # Use sequence-specific comment if available
                if comments and i - 1 < len(comments) and comments[i - 1]:
                    header += f" {comments[i-1]}"
                elif comment:  # Fall back to global comment
                    header += f" {comment}"

                # Write the sequence in chunks of 80 characters for readability
                fh.write(f"{header}\n")
                for j in range(0, len(seq), 80):
                    fh.write(f"{seq[j:j+80]}\n")
        logging.info("FASTA file written to %s", filename)
    except Exception as e:
        logging.error("Error writing FASTA file: %s", e)
        raise

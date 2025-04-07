# muc_one_up/fasta_writer.py

import logging


def write_fasta(sequences, filename, prefix="haplotype", comment=None):
    """
    Write a list of sequences to a FASTA file.

    :param sequences: List of sequence strings.
    :param filename: Output filename for the FASTA file.
    :param prefix: Header prefix for each sequence (default "haplotype").
    :param comment: Optional comment to append to the header line (default None).
    """
    try:
        with open(filename, "w") as fh:
            for i, seq in enumerate(sequences, start=1):
                # Include comment as part of the header line if provided
                header = f">{prefix}_{i}"
                if comment:
                    header += f" {comment}"
                fh.write(f"{header}\n{seq}\n")
        logging.info("FASTA file written to %s", filename)
    except Exception as e:
        logging.error("Error writing FASTA file: %s", e)
        raise

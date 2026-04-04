"""Template FASTA generator for PBSIM3 template mode.

Replicates an amplicon sequence N times into a multi-record FASTA file.
PBSIM3 --strategy templ produces one read per FASTA record, so the
number of copies directly controls the number of simulated reads.
"""

import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_template_fasta(
    amplicon_fasta: str,
    num_copies: int,
    output_path: str,
) -> str:
    """Replicate an amplicon sequence into a multi-record template FASTA.

    Args:
        amplicon_fasta: Path to input FASTA with a single amplicon sequence.
        num_copies: Number of copies to generate (one read per copy).
        output_path: Path for the output template FASTA.

    Returns:
        Path to the output template FASTA file.

    Raises:
        ValueError: If num_copies < 1 or input FASTA doesn't contain
            exactly one record.
    """
    if num_copies < 1:
        raise ValueError(f"num_copies must be >= 1, got {num_copies}")

    # Read the amplicon sequence — must contain exactly one record
    records = list(SeqIO.parse(amplicon_fasta, "fasta"))
    if len(records) != 1:
        raise ValueError(
            f"amplicon_fasta must contain exactly one FASTA record, "
            f"found {len(records)}: {amplicon_fasta}"
        )
    amplicon_seq = str(records[0].seq)

    # Create output directory if needed
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Write N copies as separate FASTA records
    template_records = [
        SeqRecord(
            Seq(amplicon_seq),
            id=f"amplicon_copy_{i:03d}",
            description="",
        )
        for i in range(1, num_copies + 1)
    ]

    SeqIO.write(template_records, str(output), "fasta")

    logging.info(
        "Generated template FASTA with %d copies of %d bp amplicon: %s",
        num_copies,
        len(amplicon_seq),
        output_path,
    )

    return str(output)

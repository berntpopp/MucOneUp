"""Primer-based amplicon extraction from haplotype FASTA files.

Extracts the PCR amplicon region from a haplotype sequence by locating
forward and reverse primer binding sites. Used by the amplicon simulation
pipeline to define the template region for PBSIM3.
"""

import logging
from pathlib import Path
from typing import NamedTuple

from ...exceptions import AmpliconExtractionError
from .sequence_utils import find_primer_binding_sites, reverse_complement


class AmpliconResult(NamedTuple):
    """Result of amplicon extraction.

    Attributes:
        sequence: Extracted amplicon DNA sequence (includes primers).
        length: Length of extracted amplicon in bp.
        forward_pos: 0-based start position of forward primer in original sequence.
        reverse_pos: 0-based start position of reverse primer RC in original sequence.
        fasta_path: Path to output FASTA file containing the amplicon.
    """

    sequence: str
    length: int
    forward_pos: int
    reverse_pos: int
    fasta_path: str


class AmpliconExtractor:
    """Extract amplicon region from a haplotype FASTA using primer sequences.

    Finds forward and reverse primer binding sites, extracts the region
    between them (inclusive of primer sequences), and writes to a new FASTA.

    Args:
        forward_primer: Forward primer sequence (5'->3').
        reverse_primer: Reverse primer sequence (5'->3', will be RC'd for search).
        expected_product_range: Optional (min, max) inclusive size range in bp.
    """

    def __init__(
        self,
        forward_primer: str,
        reverse_primer: str,
        expected_product_range: tuple[int, int] | None = None,
    ) -> None:
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.expected_product_range = expected_product_range
        self._reverse_primer_rc = reverse_complement(reverse_primer)

    def extract(self, haplotype_fasta: str, output_path: str) -> AmpliconResult:
        """Extract amplicon region from a haplotype FASTA file.

        Args:
            haplotype_fasta: Path to input haplotype FASTA (single sequence).
            output_path: Path for output FASTA containing the extracted amplicon.

        Returns:
            AmpliconResult with extracted sequence details.

        Raises:
            AmpliconExtractionError: If primers not found, ambiguous, or
                amplicon outside expected size range.
        """
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        records = list(SeqIO.parse(haplotype_fasta, "fasta"))
        if not records:
            raise AmpliconExtractionError(f"No sequences found in {haplotype_fasta}")

        record = records[0]
        template = str(record.seq)
        seq_id = record.id

        # Find forward primer binding sites
        fwd_sites = find_primer_binding_sites(template, self.forward_primer)
        if not fwd_sites:
            raise AmpliconExtractionError(
                f"Forward primer '{self.forward_primer[:20]}...' not found in "
                f"{seq_id}. Ensure the simulated reference includes sufficient "
                f"flanking sequence to contain the primer binding site."
            )
        if len(fwd_sites) > 1:
            raise AmpliconExtractionError(
                f"Multiple forward primer binding sites found at positions "
                f"{fwd_sites} in {seq_id}. Amplicon extraction requires "
                f"unambiguous primer binding."
            )

        # Find reverse primer binding sites (search for RC on forward strand)
        rev_sites = find_primer_binding_sites(
            template, self.reverse_primer, reverse_complement=True
        )
        if not rev_sites:
            raise AmpliconExtractionError(
                f"Reverse primer '{self.reverse_primer[:20]}...' not found in "
                f"{seq_id}. Ensure the simulated reference includes sufficient "
                f"flanking sequence to contain the primer binding site."
            )
        if len(rev_sites) > 1:
            raise AmpliconExtractionError(
                f"Multiple reverse primer binding sites found at positions "
                f"{rev_sites} in {seq_id}. Amplicon extraction requires "
                f"unambiguous primer binding."
            )

        fwd_pos = fwd_sites[0]
        rev_pos = rev_sites[0]

        # Extract amplicon: from forward primer start to reverse primer end (inclusive)
        amplicon_end = rev_pos + len(self._reverse_primer_rc)

        # Validate primer orientation
        if fwd_pos >= amplicon_end:
            raise AmpliconExtractionError(
                f"Forward primer at position {fwd_pos} is downstream of reverse primer "
                f"at position {rev_pos} in {seq_id}. Check primer orientation."
            )

        amplicon_seq = template[fwd_pos:amplicon_end]
        amplicon_len = len(amplicon_seq)

        # Validate size range if configured
        if self.expected_product_range is not None:
            min_size, max_size = self.expected_product_range
            if amplicon_len < min_size or amplicon_len > max_size:
                raise AmpliconExtractionError(
                    f"Extracted amplicon length ({amplicon_len} bp) is outside "
                    f"expected product range [{min_size}, {max_size}] bp."
                )

        # Write amplicon to output FASTA
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)

        amplicon_record = SeqRecord(
            Seq(amplicon_seq),
            id=f"{seq_id}_amplicon",
            description=f"amplicon fwd={fwd_pos} rev={rev_pos} len={amplicon_len}",
        )
        SeqIO.write([amplicon_record], str(output), "fasta")

        logging.info(
            "Extracted amplicon from %s: %d bp (fwd=%d, rev=%d)",
            seq_id,
            amplicon_len,
            fwd_pos,
            rev_pos,
        )

        return AmpliconResult(
            sequence=amplicon_seq,
            length=amplicon_len,
            forward_pos=fwd_pos,
            reverse_pos=rev_pos,
            fasta_path=str(output),
        )

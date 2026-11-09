#!/usr/bin/env python3
"""
Wrapper module for reseq tool operations.

This module provides wrapper functions for the reseq toolkit operations:
- replace_Ns: Replace N nucleotides in FASTA files
- generate_systematic_errors: Create systematic error profiles
- create_reads: Generate reads from fragments using seqToIllumina
"""

import gzip
import logging
import re
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command
from ..utils.tool_version import get_tool_version


def replace_Ns(
    input_fa: str,
    output_fa: str,
    tools: dict[str, str],
    seed: int | None = None,
) -> None:
    """
    Replace Ns in the simulated FASTA using reseq replaceN.

    Args:
        input_fa: Input FASTA filename.
        output_fa: Output FASTA filename.
        tools: Dictionary of tool commands.
        seed: Random seed for reproducibility.

    Raises:
        SystemExit: If the command fails.
    """
    # From reseq help: replaceN expects -r for input and -R for output
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd_args: list[str | int] = ["replaceN", "-r", input_fa, "-R", output_fa]
    if seed is not None:
        cmd_args.extend(["--seed", seed])
    cmd = build_tool_command(tools["reseq"], *cmd_args)
    run_command(cmd, timeout=60, stderr_log_level=logging.INFO, stderr_prefix="[reseq] ")


def generate_systematic_errors(
    input_fa: str,
    reseq_model: str,
    output_fq: str,
    tools: dict[str, str],
    seed: int | None = None,
) -> None:
    """
    Generate systematic errors using reseq illuminaPE.

    Args:
        input_fa: Input FASTA filename.
        reseq_model: Reseq model file.
        output_fq: Output FASTQ filename.
        tools: Dictionary of tool commands.
        seed: Random seed for reproducibility.

    Raises:
        SystemExit: If the command fails or the output file is not created.
    """
    # Based on the error and original code, this needs to be exactly like the original
    # ERROR: Either bamIn or statsIn option are mandatory.
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd_args: list[str | int] = [
        "illuminaPE",
        "-r",
        input_fa,  # Input reference fasta
        "-s",
        reseq_model,  # Reseq model file (statsIn)
        "--stopAfterEstimation",
        "--writeSysError",
        output_fq,  # Output FASTQ filename
    ]
    if seed is not None:
        cmd_args.extend(["--seed", seed])
    cmd = build_tool_command(tools["reseq"], *cmd_args)
    run_command(cmd, timeout=60, stderr_log_level=logging.INFO, stderr_prefix="[reseq] ")

    # Verify output exists and is non-empty
    output_path = Path(output_fq)
    if not output_path.exists() or output_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to generate systematic errors: Output {output_fq} missing or empty"
        )


# Legacy ReSeq 1.x (the bioconda ``reseq`` package) deadlocks in seqToIllumina on
# more than 10,000 fragments; ReSeq2 (github.com/berntpopp/ReSeq2) fixed this.
RESEQ_MAX_FRAGMENTS = 10_000


def _is_legacy_reseq(reseq_cmd: str) -> bool:
    """True only if the tool reports a 1.x version; unknown ("N/A") counts as ReSeq2 (shipped)."""
    version = get_tool_version(reseq_cmd, "reseq")
    legacy = version != "N/A" and re.search(r"version\s+([2-9]|\d{2,})", version) is None
    if version == "N/A":
        logging.warning("ReSeq version unknown; assuming ReSeq2 (no input chunking)")
    else:
        logging.info(
            "Detected %s; using %s mode", version, "legacy chunked" if legacy else "ReSeq2"
        )
    return legacy


def _split_fasta(path: str, max_records: int, prefix: str) -> list[str]:
    """Split a FASTA into files of at most *max_records* records; return their paths."""
    chunks: list[str] = []
    out = None
    count = 0
    try:
        with Path(path).open() as fin:
            for line in fin:
                if line.startswith(">"):
                    if out is None or count == max_records:
                        if out is not None:
                            out.close()
                        chunks.append(f"{prefix}.chunk{len(chunks)}.fa")
                        out = Path(chunks[-1]).open("w")
                        count = 0
                    count += 1
                if out is not None:
                    out.write(line)
    finally:
        if out is not None:
            out.close()
    return chunks


def create_reads(
    input_fragments: str,
    reseq_model: str,
    output_reads: str,
    threads: int,
    tools: dict[str, str],
    timeout: int | None = None,
    seed: int | None = None,
) -> None:
    """
    Create reads from fragments using reseq/reseq2 seqToIllumina.

    With legacy ReSeq 1.x, inputs larger than ``RESEQ_MAX_FRAGMENTS`` are
    simulated chunk by chunk (chunk *i* uses ``seed + i``) and concatenated to
    avoid its deadlock. ReSeq2 always gets the whole input in one call.

    Args:
        input_fragments: Input fragments FASTA.
        reseq_model: Reseq model file.
        output_reads: Output reads FASTQ.
        threads: Number of threads.
        tools: Dictionary of tool commands.
        timeout: Timeout in seconds per reseq call (default: none).
        seed: Random seed for reproducibility.

    Raises:
        FileOperationError: If the command fails or the output file is not created.
    """
    chunks = (
        _split_fasta(input_fragments, RESEQ_MAX_FRAGMENTS, output_reads)
        if _is_legacy_reseq(tools["reseq"])
        else []
    )
    single = len(chunks) <= 1
    if single:
        for c in chunks:
            Path(c).unlink()
        chunks = [input_fragments]
    try:
        for i, chunk in enumerate(chunks):
            chunk_out = output_reads if single else f"{chunk}.fq"
            cmd_args: list[str | int] = [
                "seqToIllumina",
                "-j",
                threads,
                "-s",
                reseq_model,
                "-i",
                chunk,
                "-o",
                chunk_out,
            ]
            if seed is not None:
                cmd_args.extend(["--seed", seed + i])
            run_command(
                build_tool_command(tools["reseq"], *cmd_args),
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[reseq] ",
            )
        if not single:
            with Path(output_reads).open("wb") as fout:
                for chunk in chunks:
                    chunk_fq = Path(f"{chunk}.fq")
                    if chunk_fq.exists():
                        fout.write(chunk_fq.read_bytes())
    finally:
        if not single:
            for chunk in chunks:
                Path(chunk).unlink(missing_ok=True)
                Path(f"{chunk}.fq").unlink(missing_ok=True)

    # Verify output exists and is non-empty
    output_path = Path(output_reads)
    if not output_path.exists() or output_path.stat().st_size == 0:
        raise FileOperationError(f"Failed to create reads: Output {output_reads} missing or empty")


def _template_name(header: str) -> str:
    """Read name without the comment and a trailing /1 or /2 mate suffix."""
    name = header[1:].split(maxsplit=1)[0] if header.startswith("@") and header.strip() else ""
    return name.removesuffix("/1").removesuffix("/2")


def split_reads(interleaved_fastq: str, output_fastq1: str, output_fastq2: str) -> None:
    """
    Split an interleaved FASTQ (4 lines per record) into two gzipped FASTQ files.

    Args:
        interleaved_fastq: Input interleaved FASTQ filename.
        output_fastq1: Output FASTQ filename for read1.
        output_fastq2: Output FASTQ filename for read2.

    Raises:
        FileOperationError: If the input is truncated, has an odd number of
            records, or an R2 header does not match its R1 header.
        OSError: If the input file can't be read or output files can't be written.
    """
    with (
        Path(interleaved_fastq).open() as in_fq,
        gzip.open(output_fastq1, "wt", compresslevel=1) as out_fq1,
        gzip.open(output_fastq2, "wt", compresslevel=1) as out_fq2,
    ):
        pair = 0
        while True:
            read1 = [in_fq.readline() for _ in range(4)]
            if not read1[0]:
                break
            pair += 1
            read2 = [in_fq.readline() for _ in range(4)]
            if not all(read1 + read2):  # readline() returns "" at EOF
                raise FileOperationError(
                    f"Truncated interleaved FASTQ {interleaved_fastq}: pair {pair} is incomplete"
                )
            name1, name2 = _template_name(read1[0]), _template_name(read2[0])
            if not name1 or name1 != name2:
                raise FileOperationError(
                    f"Interleaved FASTQ {interleaved_fastq} out of step at pair {pair}: "
                    f"R1 {read1[0].strip()!r} vs R2 {read2[0].strip()!r}"
                )
            out_fq1.writelines(read1)
            out_fq2.writelines(read2)

    logging.info(
        "Split interleaved FASTQ %s into %s and %s",
        interleaved_fastq,
        output_fastq1,
        output_fastq2,
    )

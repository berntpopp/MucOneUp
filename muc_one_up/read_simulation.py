#!/usr/bin/env python3
"""
read_simulation.py

This module integrates a read‐simulation pipeline into the muconeup project.
It uses external tools (reseq, faToTwoBit, samtools, pblat, bwa) and a ported version
of w‑Wessim2 (originally written in Python 2) to simulate Illumina reads from a simulated
MUC1 haplotype FASTA.

The pipeline performs the following steps:
 1. Replace Ns in the simulated FASTA using reseq.
 2. Generate systematic errors (illuminaPE) to create an error profile.
 3. Convert the “noNs” FASTA to 2bit format.
 4. Extract a subset reference from a sample BAM using samtools.
 5. Run pblat to align the 2bit file to the subset reference.
 6. Simulate fragments (port of w‑Wessim2) from the “noNs” FASTA.
 7. Create reads from the fragments using reseq seqToIllumina.
 8. Split the interleaved FASTQ into paired FASTQ files (gzipped).
 9. Align the reads to a human reference with bwa mem and sort/index with samtools.
10. Optionally, downsample the VNTR region in the final BAM file to a fixed target coverage.
11. Clean up all intermediate files.

Usage:
    python read_simulation.py <config.json> <input_fasta>
where <input_fasta> is typically the output from muconeup (e.g. muc1_simulated.fa).
"""

import gzip
import logging
import math
import os
import random
import signal
import subprocess
import sys
import time
from datetime import datetime
import threading
import shutil  # NEW: for checking executables in PATH
from pathlib import Path  # NEW: for downsample_bam

# Configure logging.
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


def run_command(cmd, shell=False, timeout=None):
    """
    Run a command in its own process group so that it can be killed on timeout.
    Capture both stdout and stderr live and log them line-by-line.

    :param cmd: Command (list or string) to run.
    :param shell: Whether to run the command in the shell.
    :param timeout: Timeout in seconds (None means wait indefinitely).
    :return: The process return code.
    """
    if isinstance(cmd, list) and not shell and " " in cmd[0]:
        shell = True
        cmd = " ".join(cmd)
    cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
    logging.info("Running command: %s (timeout=%s)", cmd_str, timeout)

    try:
        proc = subprocess.Popen(
            cmd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            preexec_fn=os.setsid
        )
    except Exception as e:
        logging.exception("Failed to start command: %s", cmd_str)
        sys.exit(1)

    def log_stream(stream, log_func):
        """Read lines from a stream and log them using log_func."""
        for line in iter(stream.readline, ""):
            log_func(line.rstrip())
        stream.close()

    stdout_thread = threading.Thread(target=log_stream, args=(proc.stdout, logging.info))
    stderr_thread = threading.Thread(target=log_stream, args=(proc.stderr, logging.error))
    stdout_thread.start()
    stderr_thread.start()

    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        logging.warning("Command timed out after %s seconds. Killing process group.", timeout)
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except Exception as e:
            logging.exception("Error killing process group")
        proc.wait()
    stdout_thread.join()
    stderr_thread.join()

    if proc.returncode != 0:
        logging.error("Command exited with non-zero exit code %d: %s", proc.returncode, cmd_str)
        sys.exit(proc.returncode)
    return proc.returncode


def fix_field(field, desired_length, pad_char):
    """
    Ensure a field is exactly desired_length characters long.

    :param field: Input string.
    :param desired_length: Desired length.
    :param pad_char: Character to pad with.
    :return: String of length desired_length.
    """
    if len(field) > desired_length:
        return field[:desired_length]
    else:
        return field.ljust(desired_length, pad_char)


def replace_Ns(input_fa, output_fa, tools):
    """
    Replace Ns in the simulated FASTA using reseq replaceN.

    :param input_fa: Input FASTA filename.
    :param output_fa: Output FASTA filename.
    :param tools: Dictionary of tool commands.
    """
    cmd = [tools["reseq"], "replaceN", "-r", input_fa, "-R", output_fa]
    run_command(cmd, timeout=60)


def generate_systematic_errors(input_fa, reseq_model, output_fq, tools):
    """
    Generate systematic errors using reseq illuminaPE.

    :param input_fa: Input FASTA filename.
    :param reseq_model: Reseq model file.
    :param output_fq: Output FASTQ filename.
    :param tools: Dictionary of tool commands.
    """
    cmd = [
        tools["reseq"], "illuminaPE", "-r", input_fa, "-s", reseq_model,
        "--stopAfterEstimation", "--writeSysError", output_fq
    ]
    run_command(cmd, timeout=60)


def fa_to_twobit(input_fa, output_2bit, tools):
    """
    Convert a FASTA file to 2bit format using faToTwoBit.

    :param input_fa: Input FASTA filename.
    :param output_2bit: Output 2bit filename.
    :param tools: Dictionary of tool commands.
    """
    cmd = [tools["faToTwoBit"], input_fa, output_2bit]
    run_command(cmd, timeout=60)


def extract_subset_reference(sample_bam, output_fa, tools):
    """
    Extract a subset reference from a BAM file.

    :param sample_bam: Input BAM filename.
    :param output_fa: Output FASTA filename.
    :param tools: Dictionary of tool commands.
    :return: Intermediate collated BAM filename.
    """
    intermediate_bam = output_fa + ".collated.bam"
    cmd_collate = f"{tools['samtools']} collate -u {sample_bam} -o {intermediate_bam}"
    run_command(cmd_collate, shell=True, timeout=60)
    cmd_fasta = f"{tools['samtools']} fasta {intermediate_bam} > {output_fa}"
    run_command(cmd_fasta, shell=True, timeout=60)
    return intermediate_bam


def run_pblat(twobit_file, subset_reference, output_psl, tools, threads=24,
              minScore=95, minIdentity=95):
    """
    Run pblat to align a 2bit file to a subset reference.

    :param twobit_file: Input 2bit filename.
    :param subset_reference: Subset reference FASTA.
    :param output_psl: Output PSL filename.
    :param tools: Dictionary of tool commands.
    :param threads: Number of threads.
    :param minScore: Minimal score.
    :param minIdentity: Minimal identity.
    """
    cmd = [
        tools["pblat"], twobit_file, subset_reference, output_psl,
        f"-threads={threads}", f"-minScore={minScore}", f"-minIdentity={minIdentity}"
    ]
    run_command(cmd, timeout=120)


def read_fasta_to_dict(fasta_file):
    """
    Read a FASTA file into a dictionary mapping chromosome names to sequences.

    :param fasta_file: Path to the FASTA file.
    :return: Dictionary mapping chromosome names to sequences.
    """
    ref_dict = {}
    chrom = None
    seq_lines = []
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if chrom is not None:
                    ref_dict[chrom] = "".join(seq_lines)
                chrom = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if chrom is not None:
            ref_dict[chrom] = "".join(seq_lines)
    return ref_dict


def read_syser_file(syser_file):
    """
    Read the systematic errors file and return four dictionaries.

    :param syser_file: Path to the syser file.
    :return: Tuple (tendendict_f, tendendict_r, ratesdict_f, ratesdict_r).
    """
    tendendict_f = {}
    tendendict_r = {}
    ratesdict_f = {}
    ratesdict_r = {}
    direction = None
    rat = False
    with open(syser_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("@") and len(line) < 500:
                space_index = line.find(" ")
                identifier = line[1:space_index] if space_index != -1 else line[1:]
                if "forward" in line:
                    direction = "f"
                elif "reverse" in line:
                    direction = "r"
                else:
                    logging.error("Direction not found in syser header")
                continue
            elif line == "+":
                rat = True
                continue
            else:
                if rat:
                    if direction == "f":
                        ratesdict_f[identifier] = line
                    elif direction == "r":
                        ratesdict_r[identifier] = line[::-1]
                    rat = False
                else:
                    if direction == "f":
                        tendendict_f[identifier] = line
                    elif direction == "r":
                        tendendict_r[identifier] = line[::-1]
    return tendendict_f, tendendict_r, ratesdict_f, ratesdict_r


def read_psl_file(psl_file):
    """
    Read a PSL file (skipping header lines) and return a list of matches.

    Each match is a tuple: (chrom, start, end, strand).

    :param psl_file: Path to the PSL file.
    :return: List of match tuples.
    """
    matches = []
    with open(psl_file, "r") as fh:
        for _ in range(5):
            fh.readline()
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 17:
                continue
            try:
                qgapsize = int(parts[5])
                tgapsize = int(parts[7])
            except ValueError:
                continue
            if qgapsize > 2 or tgapsize > 2:
                continue
            chrom = parts[13]
            try:
                start = int(parts[15])
                end = int(parts[16])
            except ValueError:
                continue
            strand = parts[8]
            matches.append((chrom, start, end, strand))
    logging.debug("Found %d matches in PSL file.", len(matches))
    return matches


def pick_on_match(matches):
    """
    Randomly pick one match from the provided list.

    :param matches: List of match tuples.
    :return: A single match tuple.
    """
    match = random.choice(matches)
    logging.debug("Picked match: %s", match)
    return match


def get_insert_length(mu, sigma, lower):
    """
    Return an insert length sampled from a Gaussian distribution (>= lower).

    :param mu: Mean insert length.
    :param sigma: Standard deviation.
    :param lower: Minimum allowed length.
    :return: Integer insert length.
    """
    while True:
        length = int(random.gauss(mu, sigma))
        if length >= lower:
            logging.debug("Selected insert length: %d", length)
            return length


def pick_fragment(match, ins, bind):
    """
    Randomly pick a fragment based on a match and desired insert length.

    :param match: Tuple (chrom, probe_start, probe_end, strand).
    :param ins: Desired insert length.
    :param bind: Minimum fraction (%) for overlap.
    :return: Tuple (chrom, fragment_start, fragment_end, strand).
    """
    probechrom, probestart, probeend, strand = match
    probelength = probeend - probestart
    minimummatch = int(probelength * bind / 100)
    overlap = int(random.triangular(minimummatch, probelength, probelength))
    margin = max(ins - overlap, 0)
    rangestart = probestart - margin
    rangeend = probeend + margin
    if rangeend - ins < rangestart:
        seqstart = rangestart
    else:
        seqstart = random.randint(rangestart, rangeend - ins)
    fragment = (probechrom, seqstart, seqstart + ins, strand)
    logging.debug("Picked fragment: %s", fragment)
    return fragment


def comp(sequence):
    """
    Return the complement of a DNA sequence (preserving case).

    :param sequence: DNA sequence.
    :return: Complemented sequence.
    """
    d = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "a": "t", "t": "a", "c": "g", "g": "c",
        "N": "N", "n": "n"
    }
    return "".join(d.get(s, "N") for s in sequence)


# NEW: External tool check function
def check_external_tools(tools):
    """
    Check that all external tools required for read simulation are available.

    For each tool (e.g. reseq, faToTwoBit, samtools, pblat, bwa), this function:
      - Validates that the executable is in the system PATH.
      - Runs a basic test command (using '--help' or '--version') for tools that support it.
        Currently, only "reseq" and "samtools" are tested.
      - For tools that do not support a test parameter (e.g. bwa, pblat, faToTwoBit), only verifies they are found.
    
    Exits gracefully with a clear error message if any required tool is missing or misconfigured.
    """
    test_args = {
        "reseq": "--help",
        "samtools": "--version"
    }
    
    for tool_key, cmd_str in tools.items():
        tokens = cmd_str.split()
        executable = tokens[0]
        if shutil.which(executable) is None:
            logging.error(
                f"Required tool '{tool_key}' not found: {executable}. Please ensure it is installed and in your PATH."
            )
            sys.exit(1)
        
        if tool_key in test_args:
            test_command = f"{cmd_str} {test_args[tool_key]}"
            try:
                result = subprocess.run(
                    test_command,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=30
                )
                if result.returncode != 0:
                    logging.error(
                        f"Test command for tool '{tool_key}' failed:\n"
                        f"  Command: {test_command}\n"
                        f"  Error Output: {result.stderr.decode().strip()}"
                    )
                    sys.exit(1)
                else:
                    logging.info(f"Tool '{tool_key}' is available and working.")
            except Exception as e:
                logging.error(
                    f"Exception when testing tool '{tool_key}' with command '{test_command}': {e}"
                )
                sys.exit(1)
        else:
            logging.info(f"Tool '{tool_key}' found: {executable}")


# NEW: Function to calculate mean VNTR coverage using samtools depth
def calculate_vntr_coverage(samtools_exe, bam_file, region, threads, output_dir, output_name):
    """
    Calculate mean coverage over the VNTR region using samtools depth.

    :param samtools_exe: Path to the samtools executable.
    :param bam_file: BAM file for which coverage is calculated.
    :param region: Genomic region in 'chr:start-end' format.
    :param threads: Number of threads.
    :param output_dir: Directory to write the coverage output.
    :param output_name: Base name for the coverage output file.
    :return: Mean coverage (float).
    """
    coverage_output = os.path.join(output_dir, f"{output_name}_vntr_coverage.txt")
    depth_command = [
        samtools_exe,
        "depth",
        "-@",
        str(threads),
        "-r",
        region,
        str(bam_file)
    ]
    logging.info("Calculating VNTR coverage for %s in region %s", bam_file, region)
    with open(coverage_output, "w") as fout:
        if " " in samtools_exe:
            cmd_str = " ".join(depth_command)
            subprocess.run(cmd_str, shell=True, stdout=fout, check=True)
        else:
            subprocess.run(depth_command, stdout=fout, check=True)
    coverage_values = []
    with open(coverage_output, "r") as fin:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                try:
                    coverage_values.append(int(parts[2]))
                except ValueError:
                    continue
    if coverage_values:
        mean_cov = sum(coverage_values) / len(coverage_values)
    else:
        mean_cov = 0
    logging.info("Mean VNTR coverage: %.2f", mean_cov)
    return mean_cov


# NEW: Function to calculate mean target (non-VNTR) coverage using samtools depth with BED file
def calculate_target_coverage(samtools_exe, bam_file, bed_file, threads, output_dir, output_name):
    """
    Calculate mean coverage over the target (non-VNTR) regions using samtools depth and a BED file.

    :param samtools_exe: Path to the samtools executable.
    :param bam_file: BAM file for which coverage is calculated.
    :param bed_file: BED file specifying the target regions.
    :param threads: Number of threads.
    :param output_dir: Directory to write the coverage output.
    :param output_name: Base name for the coverage output file.
    :return: Mean coverage (float).
    """
    coverage_output = os.path.join(output_dir, f"{output_name}_target_coverage.txt")
    depth_command = [
        samtools_exe,
        "depth",
        "-@",
        str(threads),
        "-b",
        bed_file,
        str(bam_file)
    ]
    logging.info("Calculating target (non-VNTR) coverage for %s using bed file %s", bam_file, bed_file)
    with open(coverage_output, "w") as fout:
        if " " in samtools_exe:
            cmd_str = " ".join(depth_command)
            subprocess.run(cmd_str, shell=True, stdout=fout, check=True)
        else:
            subprocess.run(depth_command, stdout=fout, check=True)
    coverage_values = []
    with open(coverage_output, "r") as fin:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                try:
                    coverage_values.append(int(parts[2]))
                except ValueError:
                    continue
    if coverage_values:
        mean_cov = sum(coverage_values) / len(coverage_values)
    else:
        mean_cov = 0
    logging.info("Mean target coverage: %.2f", mean_cov)
    return mean_cov


# NEW: Function to downsample BAM file using samtools view -s for a specified region
def downsample_bam(samtools_exe, input_bam, output_bam, region, fraction, seed, threads):
    """
    Downsample the BAM file to the specified fraction in the given region.
    """
    output_bam = Path(output_bam)
    subsample_param = f"{seed}.{int(fraction * 1000):03d}"
    view_command = [
        samtools_exe,
        "view",
        "-s",
        subsample_param,
        "-@",
        str(threads),
        "-b",
        "-o",
        str(output_bam),
        str(input_bam),
        region,
    ]
    run_command(view_command, timeout=60)
    sorted_bam = output_bam.with_suffix(".sorted.bam")
    sort_command = [
        samtools_exe,
        "sort",
        "-@",
        str(threads),
        "-o",
        str(sorted_bam),
        str(output_bam)
    ]
    run_command(sort_command, timeout=60)
    os.remove(str(output_bam))
    os.rename(str(sorted_bam), str(output_bam))
    index_command = [samtools_exe, "index", str(output_bam)]
    run_command(index_command, timeout=60)


# NEW: Function to downsample entire BAM file (without region restriction) using samtools view -s
def downsample_entire_bam(samtools_exe, input_bam, output_bam, fraction, seed, threads):
    """
    Downsample the entire BAM file to the specified fraction.
    """
    output_bam = Path(output_bam)
    subsample_param = f"{seed}.{int(fraction * 1000):03d}"
    view_command = [
        samtools_exe,
        "view",
        "-s",
        subsample_param,
        "-@",
        str(threads),
        "-b",
        "-o",
        str(output_bam),
        str(input_bam)
    ]
    run_command(view_command, timeout=60)
    sorted_bam = output_bam.with_suffix(".sorted.bam")
    sort_command = [
        samtools_exe,
        "sort",
        "-@",
        str(threads),
        "-o",
        str(sorted_bam),
        str(output_bam)
    ]
    run_command(sort_command, timeout=60)
    os.remove(str(output_bam))
    os.rename(str(sorted_bam), str(output_bam))
    index_command = [samtools_exe, "index", str(output_bam)]
    run_command(index_command, timeout=60)


def simulate_fragments(ref_fa, syser_file, psl_file, read_number, fragment_size,
                       fragment_sd, min_fragment, bind, output_fragments):
    """
    Simulate fragments (port of w‑Wessim2) and write paired fragment sequences to a FASTA file.

    :param ref_fa: Reference FASTA file.
    :param syser_file: Systematic errors file.
    :param psl_file: PSL file.
    :param read_number: Number of fragments to simulate.
    :param fragment_size: Desired fragment size.
    :param fragment_sd: Standard deviation for fragment size.
    :param min_fragment: Minimum acceptable fragment length.
    :param bind: Minimum fraction (%) for overlap.
    :param output_fragments: Output FASTA filename for fragments.
    """
    t0 = time.time()
    logging.info("Simulating fragments...")
    ref_dict = read_fasta_to_dict(ref_fa)
    tendendict_f, tendendict_r, ratesdict_f, ratesdict_r = read_syser_file(syser_file)
    matches = read_psl_file(psl_file)
    if not matches:
        logging.error("No matches found in PSL file.")
        sys.exit(1)
    output_dir = os.path.dirname(output_fragments)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_fragments, "w") as out_f:
        i = 0
        attempts = 0
        max_attempts = read_number * 10
        while i < read_number and attempts < max_attempts:
            attempts += 1
            match = pick_on_match(matches)
            ins = get_insert_length(fragment_size, fragment_sd, min_fragment)
            fragment = pick_fragment(match, ins, bind)
            chrom, frag_start, frag_end, strand = fragment
            seq = ref_dict.get(chrom, "")[frag_start:frag_end]
            if len(seq) < min_fragment:
                logging.debug("Fragment length %d is below min_fragment %d. Skipping.",
                              len(seq), min_fragment)
                continue
            frag_seq = seq
            frag_seq_rc = comp(seq)[::-1]
            frag_len = len(frag_seq)
            tendenseq_f = tendendict_f.get(chrom, "N" * frag_len)
            ratesseq_f = ratesdict_f.get(chrom, "!" * frag_len)
            tendenseq_r = tendendict_r.get(chrom, "N" * frag_len)
            ratesseq_r = ratesdict_r.get(chrom, "!" * frag_len)
            tendenseq_f_fixed = fix_field(tendenseq_f, frag_len, "N")
            ratesseq_f_fixed = fix_field(ratesseq_f, frag_len, "!")
            tendenseq_r_fixed = fix_field(tendenseq_r, frag_len, "N")
            ratesseq_r_fixed = fix_field(ratesseq_r, frag_len, "!")
            if strand == "+":
                out_f.write(
                    f">{i+1} 1;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n{frag_seq}\n"
                )
                out_f.write(
                    f">{i+1} 2;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n{frag_seq_rc}\n"
                )
            elif strand == "-":
                out_f.write(
                    f">{i+1} 1;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n{frag_seq_rc}\n"
                )
                out_f.write(
                    f">{i+1} 2;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n{frag_seq}\n"
                )
            i += 1
            if (i + 1) % 1000000 == 0:
                t1 = time.time()
                logging.info("%d fragments processed in %.2f secs", i + 1, t1 - t0)
        if i < read_number:
            logging.warning("Only generated %d fragments after %d attempts. Expected %d.",
                            i, attempts, read_number)
    t1 = time.time()
    logging.info("Done processing %d fragments in %.2f secs", i, t1 - t0)


def create_reads(input_fragments, reseq_model, output_reads, threads, tools):
    """
    Create reads from fragments using reseq seqToIllumina.

    If the command times out but the output file exists and is non-empty,
    a warning is logged and the process continues.

    :param input_fragments: Input fragments FASTA.
    :param reseq_model: Reseq model file.
    :param output_reads: Output reads FASTQ.
    :param threads: Number of threads.
    :param tools: Dictionary of tool commands.
    """
    cmd = [
        tools["reseq"], "seqToIllumina", "-j", str(threads),
        "-s", reseq_model, "-i", input_fragments, "-o", output_reads
    ]
    try:
        run_command(cmd, timeout=120)
    except SystemExit as e:
        if os.path.exists(output_reads) and os.path.getsize(output_reads) > 0:
            logging.warning(
                "seqToTwoIllumina command timed out but output file exists. Continuing."
            )
        else:
            raise


def split_reads(interleaved_fastq, output_fastq1, output_fastq2):
    """
    Split an interleaved FASTQ (4 lines per record) into two gzipped FASTQ files.

    :param interleaved_fastq: Input interleaved FASTQ filename.
    :param output_fastq1: Output FASTQ filename for read1.
    :param output_fastq2: Output FASTQ filename for read2.
    """
    logging.info("Splitting %s into %s and %s", interleaved_fastq, output_fastq1, output_fastq2)
    with open(interleaved_fastq, "r") as inf, \
         gzip.open(output_fastq1, "wt") as out1, \
         gzip.open(output_fastq2, "wt") as out2:
        record = []
        record_index = 0
        for line in inf:
            record.append(line)
            if len(record) == 4:
                if record_index % 2 == 0:
                    out1.writelines(record)
                else:
                    out2.writelines(record)
                record = []
                record_index += 1
    logging.info("Splitting complete.")


def align_reads(read1, read2, human_reference, output_bam, tools, threads=4):
    """
    Align paired-end reads with bwa mem, sort with samtools, and index the BAM file.

    :param read1: FASTQ filename for read1.
    :param read2: FASTQ filename for read2.
    :param human_reference: Human reference FASTA.
    :param output_bam: Output BAM filename.
    :param tools: Dictionary of tool commands.
    :param threads: Number of threads.
    """
    bwa_exe = tools["bwa"]
    if " " in bwa_exe:
        bwa_cmd = "bash -c '" + " ".join([bwa_exe, "mem", human_reference, read1, read2]) + "'"
        bwa_shell = True
    else:
        bwa_cmd = [bwa_exe, "mem", human_reference, read1, read2]
        bwa_shell = False

    samtools_exe = tools["samtools"]
    if " " in samtools_exe:
        samtools_sort_cmd = (
            "bash -c '" +
            " ".join([samtools_exe, "sort", "-o", output_bam, "-"]) + "'"
        )
        samtools_shell = True
    else:
        samtools_sort_cmd = [samtools_exe, "sort", "-o", output_bam, "-"]
        samtools_shell = False

    logging.info("Aligning reads...")
    try:
        bwa_proc = subprocess.Popen(
            bwa_cmd,
            shell=bwa_shell,
            stdout=subprocess.PIPE,
            universal_newlines=True
        )
        sort_proc = subprocess.Popen(
            samtools_sort_cmd,
            shell=samtools_shell,
            stdin=bwa_proc.stdout,
            stdout=subprocess.PIPE,
            universal_newlines=True
        )
        bwa_proc.stdout.close()
        sort_proc.communicate()
    except Exception as e:
        logging.exception("Error during alignment")
        sys.exit(1)
    index_cmd = [tools["samtools"], "index", output_bam]
    run_command(index_cmd, timeout=60)


def cleanup_files(file_list):
    """
    Remove files in the provided list if they exist.

    :param file_list: List of filenames.
    """
    for f in file_list:
        try:
            if os.path.exists(f):
                os.remove(f)
                logging.info("Removed intermediate file: %s", f)
        except Exception as e:
            logging.warning("Unable to remove intermediate file %s: %s", f, e)


def simulate_reads(config, input_fa):
    """
    Run the complete read simulation pipeline.

    :param config: Dictionary containing "tools" and "read_simulation" sections.
    :param input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).
    """
    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})

    # Check that all external tools are available and working.
    check_external_tools(tools)

    base = os.path.splitext(input_fa)[0]
    noNs_fa = base + "_noNs.fasta"
    syser_fq = base + "_noNs_syserrors.fq"
    twobit_file = noNs_fa + ".2bit"
    subset_ref = base + "_subset.fasta"
    psl_file = base + "_subset.psl"
    fragments_fa = base + "_w_wessim2.fa"
    reads_fq = base + "_w_wessim2.fq"
    reads_fq1 = base + "_w_wessim2_1.fq.gz"
    reads_fq2 = base + "_w_wessim2_2.fq.gz"
    output_bam = base + "_wessim.bam"

    reseq_model = rs_config.get("reseq_model")
    sample_bam = rs_config.get("sample_bam")
    human_reference = rs_config.get("human_reference")
    read_number = int(rs_config.get("read_number", 10000))
    fragment_size = int(rs_config.get("fragment_size", 170))
    fragment_sd = int(rs_config.get("fragment_sd", 35))
    min_fragment = int(rs_config.get("min_fragment", 20))
    threads = int(rs_config.get("threads", 24))
    bind = int(rs_config.get("bind", 50))

    logging.info("Starting read simulation pipeline at %s",
                 datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    replace_Ns(input_fa, noNs_fa, tools)
    generate_systematic_errors(noNs_fa, reseq_model, syser_fq, tools)
    fa_to_twobit(noNs_fa, twobit_file, tools)
    collated_bam = extract_subset_reference(sample_bam, subset_ref, tools)
    run_pblat(twobit_file, subset_ref, psl_file, tools, threads=threads)
    simulate_fragments(noNs_fa, syser_fq, psl_file, read_number,
                       fragment_size, fragment_sd, min_fragment, bind, fragments_fa)
    if not os.path.exists(fragments_fa) or os.path.getsize(fragments_fa) == 0:
        logging.error("Fragments file %s was not created or is empty.", fragments_fa)
        sys.exit(1)
    create_reads(fragments_fa, reseq_model, reads_fq, threads, tools)
    split_reads(reads_fq, reads_fq1, reads_fq2)
    align_reads(reads_fq1, reads_fq2, human_reference, output_bam, tools, threads=threads)

    # NEW: Optional downsampling of the final BAM file to a target coverage.
    downsample_target = rs_config.get("downsample_coverage")
    if downsample_target is not None:
        try:
            downsample_target = int(downsample_target)
        except ValueError:
            logging.error("Invalid downsample_coverage value in config.")
            sys.exit(1)
        # Normalize mode to lower-case and strip whitespace
        mode = rs_config.get("downsample_mode", "vntr").strip().lower()
        if mode == "vntr":
            reference_assembly = rs_config.get("reference_assembly", "hg19")
            vntr_region_key = f"vntr_region_{reference_assembly}"
            vntr_region = rs_config.get(vntr_region_key)
            if not vntr_region:
                logging.error(f"VNTR region not specified in config for {reference_assembly}.")
                sys.exit(1)
            current_cov = calculate_vntr_coverage(tools["samtools"], output_bam, vntr_region, threads, os.path.dirname(output_bam), os.path.basename(output_bam).replace(".bam", ""))
            region_info = vntr_region
        elif mode == "non_vntr":
            bed_file = rs_config.get("sample_target_bed")
            if not bed_file:
                logging.error("For non-VNTR downsampling, 'sample_target_bed' must be provided in config.")
                sys.exit(1)
            current_cov = calculate_target_coverage(tools["samtools"], output_bam, bed_file, threads, os.path.dirname(output_bam), os.path.basename(output_bam).replace(".bam", ""))
            region_info = f"BED file: {bed_file}"
        else:
            logging.error("Invalid downsample_mode in config; use 'vntr' or 'non_vntr'.")
            sys.exit(1)
        if current_cov > downsample_target:
            fraction = downsample_target / current_cov
            fraction = min(max(fraction, 0.0), 1.0)
            logging.info("Downsampling BAM from %.2fx to target %dx (fraction: %.4f) based on %s",
                         current_cov, downsample_target, fraction, region_info)
            downsampled_bam = os.path.join(os.path.dirname(output_bam), os.path.basename(output_bam).replace(".bam", "_downsampled.bam"))
            if mode == "vntr":
                downsample_bam(tools["samtools"], output_bam, downsampled_bam, vntr_region, fraction, rs_config.get("downsample_seed", 42), threads)
            else:
                downsample_entire_bam(tools["samtools"], output_bam, downsampled_bam, fraction, rs_config.get("downsample_seed", 42), threads)
            output_bam = downsampled_bam
        else:
            logging.info("Current coverage (%.2fx) is below the target; no downsampling performed.", current_cov)

    logging.info("Read simulation pipeline completed at %s",
                 datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("Final outputs:")
    logging.info("  Aligned and indexed BAM: %s", output_bam)
    logging.info("  Paired FASTQ files (gzipped): %s and %s", reads_fq1, reads_fq2)

    intermediates = [noNs_fa, syser_fq, twobit_file, subset_ref, psl_file, fragments_fa,
                     reads_fq, collated_bam]
    cleanup_files(intermediates)


if __name__ == "__main__":
    import json
    if len(sys.argv) != 3:
        print("Usage: python read_simulation.py <config.json> <input_fasta>")
        sys.exit(1)
    config_file = sys.argv[1]
    input_fa = sys.argv[2]
    with open(config_file, "r") as fh:
        config = json.load(fh)
    simulate_reads(config, input_fa)

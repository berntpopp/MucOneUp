#!/usr/bin/env python3
"""
read_simulation.py

This module integrates a read‐simulation pipeline into the muconeup project.
It uses external tools (reseq, faToTwoBit, samtools, pblat, bwa) and a ported version
of w‐Wessim2 (originally written in Python 2) to simulate Illumina reads from a simulated
MUC1 haplotype FASTA.

The pipeline performs the following steps:
 1. Replace Ns in the simulated FASTA using reseq.
 2. Generate systematic errors (illuminaPE) to create an error profile.
 3. Convert the “noNs” FASTA to 2bit format.
 4. Extract a subset reference from a sample BAM using samtools.
 5. Run pblat to align the 2bit file to the subset reference.
 6. Simulate fragments (port of w‑Wessim2) from the “noNs” FASTA.
 7. Create reads from the fragments using reseq seqToIllumina.
 8. Split the interleaved FASTQ into paired FASTQ files (written as gzipped files).
 9. Align the reads to a human reference with bwa mem and sort/index with samtools.
10. Clean up all intermediate files.

Usage:
    python read_simulation.py <config.json> <input_fasta>
where <input_fasta> is typically the output from muconeup (e.g. muc1_simulated.fa).
"""

import subprocess
import os
import sys
import random
import math
import time
from datetime import datetime
import gzip
import logging

# Configure logging.
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(message)s")

### Helper function to run external commands ###
def run_command(cmd, shell=False):
    """
    Run a command using subprocess.run and exit if an error occurs.
    If the first element of the command (if cmd is a list) contains spaces,
    we assume it's a compound command (e.g. "mamba run -n wessim reseq")
    and switch to shell mode.
    """
    if isinstance(cmd, list) and not shell and " " in cmd[0]:
        shell = True
        cmd = " ".join(cmd)
    cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
    logging.info(f"Running command: {cmd_str}")
    try:
        result = subprocess.run(cmd, shell=shell, check=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)
        if result.stdout:
            logging.debug(f"stdout: {result.stdout}")
        if result.stderr:
            logging.debug(f"stderr: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logging.exception(f"Error running command: {cmd_str}")
        sys.exit(1)

### Helper: Ensure a field is exactly desired_length in characters.
def fix_field(field, desired_length, pad_char):
    """
    If the provided field is longer than desired_length, trim it;
    if shorter, pad on the right with pad_char.
    """
    if len(field) > desired_length:
        return field[:desired_length]
    else:
        return field.ljust(desired_length, pad_char)

### Step 1: Replace Ns using reseq replaceN ###
def replace_Ns(input_fa, output_fa, tools):
    cmd = [tools["reseq"], "replaceN", "-r", input_fa, "-R", output_fa]
    run_command(cmd)

### Step 2: Generate systematic errors using reseq illuminaPE ###
def generate_systematic_errors(input_fa, reseq_model, output_fq, tools):
    cmd = [tools["reseq"], "illuminaPE", "-r", input_fa, "-s", reseq_model,
           "--stopAfterEstimation", "--writeSysError", output_fq]
    run_command(cmd)

### Step 3: Convert FASTA to 2bit using faToTwoBit ###
def fa_to_twobit(input_fa, output_2bit, tools):
    cmd = [tools["faToTwoBit"], input_fa, output_2bit]
    run_command(cmd)

### Step 4: Extract a subset reference from a BAM ###
def extract_subset_reference(sample_bam, output_fa, tools):
    """
    First, run samtools collate to generate an intermediate BAM file.
    Then, run samtools fasta on that file to generate the subset FASTA.
    Returns the name of the intermediate BAM file.
    """
    intermediate_bam = output_fa + ".collated.bam"
    cmd_collate = f"{tools['samtools']} collate -u {sample_bam} -o {intermediate_bam}"
    run_command(cmd_collate, shell=True)
    cmd_fasta = f"{tools['samtools']} fasta {intermediate_bam} > {output_fa}"
    run_command(cmd_fasta, shell=True)
    return intermediate_bam

### Step 5: Run pblat ###
def run_pblat(twobit_file, subset_reference, output_psl, tools, threads=24, minScore=95, minIdentity=95):
    cmd = [tools["pblat"], twobit_file, subset_reference, output_psl,
           f"-threads={threads}", f"-minScore={minScore}", f"-minIdentity={minIdentity}"]
    run_command(cmd)

### --- Functions for the w-wessim2 port (Step 6) --- ###
def read_fasta_to_dict(fasta_file):
    """Read a FASTA file into a dictionary {chrom: sequence}."""
    ref_dict = {}
    chrom = None
    seq_lines = []
    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if chrom is not None:
                    ref_dict[chrom] = ''.join(seq_lines)
                chrom = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if chrom is not None:
            ref_dict[chrom] = ''.join(seq_lines)
    return ref_dict

def read_syser_file(syser_file):
    """
    Reads the systematic errors file.
    Returns dictionaries: tendendict_f, tendendict_r, ratesdict_f, ratesdict_r.
    """
    tendendict_f = {}
    tendendict_r = {}
    ratesdict_f = {}
    ratesdict_r = {}
    direction = None
    rat = False
    with open(syser_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('@') and len(line) < 500:
                space_index = line.find(' ')
                id = line[1:space_index] if space_index != -1 else line[1:]
                if 'forward' in line:
                    direction = 'f'
                elif 'reverse' in line:
                    direction = 'r'
                else:
                    logging.error("Direction not found in syser header")
                continue
            elif line == '+':
                rat = True
                continue
            else:
                if rat:
                    if direction == 'f':
                        ratesdict_f[id] = line
                    elif direction == 'r':
                        ratesdict_r[id] = line[::-1]
                    rat = False
                else:
                    if direction == 'f':
                        tendendict_f[id] = line
                    elif direction == 'r':
                        tendendict_r[id] = line[::-1]
    return tendendict_f, tendendict_r, ratesdict_f, ratesdict_r

def read_psl_file(psl_file):
    """
    Read a PSL file (skipping header) and return a list of matches.
    Each match is a tuple: (chrom, start, end, strand)
    """
    matches = []
    with open(psl_file, 'r') as fh:
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
    logging.debug(f"Found {len(matches)} matches in PSL file.")
    return matches

def pick_on_match(matches):
    """Randomly pick one match from the list."""
    match = random.choice(matches)
    logging.debug(f"Picked match: {match}")
    return match

def get_insert_length(mu, sigma, lower):
    """Return an insert length sampled from a Gaussian (>= lower)."""
    while True:
        length = int(random.gauss(mu, sigma))
        if length >= lower:
            logging.debug(f"Selected insert length: {length}")
            return length

def pick_fragment(match, ins, bind):
    """
    Given a match (chrom, probe_start, probe_end, strand) and desired insert length,
    randomly pick a fragment.
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
    logging.debug(f"Picked fragment: {fragment}")
    return fragment

def comp(sequence):
    """Return the complement of a DNA sequence (preserving case)."""
    d = {'A':'T','T':'A','C':'G','G':'C',
         'a':'t','t':'a','c':'g','g':'c',
         'N':'N','n':'n'}
    return ''.join(d.get(s, 'N') for s in sequence)

def simulate_fragments(ref_fa, syser_file, psl_file, read_number, fragment_size,
                       fragment_sd, min_fragment, bind, output_fragments):
    """
    Simulate fragments (port of w‑Wessim2).
    Reads the reference FASTA, syser file, and BLAT PSL file; then writes simulated
    paired fragment sequences in FASTA format to output_fragments.
    The header for each fragment includes two fields (systematic errors and quality)
    that must be exactly the same length as the fragment sequence.
    """
    t0 = time.time()
    logging.info("Simulating fragments...")
    ref_dict = read_fasta_to_dict(ref_fa)
    tendendict_f, tendendict_r, ratesdict_f, ratesdict_r = read_syser_file(syser_file)
    matches = read_psl_file(psl_file)
    if not matches:
        logging.error("No matches found in PSL file.")
        sys.exit(1)
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_fragments)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_fragments, 'w') as out_f:
        i = 0
        attempts = 0
        max_attempts = read_number * 10  # Avoid infinite looping
        while i < read_number and attempts < max_attempts:
            attempts += 1
            match = pick_on_match(matches)
            ins = get_insert_length(fragment_size, fragment_sd, min_fragment)
            fragment = pick_fragment(match, ins, bind)
            chrom, frag_start, frag_end, strand = fragment
            seq = ref_dict.get(chrom, "")[frag_start:frag_end]
            if len(seq) < min_fragment:
                logging.debug(f"Fragment length {len(seq)} is below min_fragment {min_fragment}. Skipping.")
                continue
            # Use the full fragment sequence
            frag_seq = seq
            frag_seq_rc = comp(seq)[::-1]
            frag_len = len(frag_seq)
            # Get systematic error strings from the dictionaries.
            tendenseq_f = tendendict_f.get(chrom, "N" * frag_len)
            ratesseq_f = ratesdict_f.get(chrom, "!" * frag_len)
            tendenseq_r = tendendict_r.get(chrom, "N" * frag_len)
            ratesseq_r = ratesdict_r.get(chrom, "!" * frag_len)
            # Ensure they are exactly the same length as the fragment.
            tendenseq_f_fixed = fix_field(tendenseq_f, frag_len, "N")
            ratesseq_f_fixed = fix_field(ratesseq_f, frag_len, "!")
            tendenseq_r_fixed = fix_field(tendenseq_r, frag_len, "N")
            ratesseq_r_fixed = fix_field(ratesseq_r, frag_len, "!")
            # Write header and sequence.
            if strand == '+':
                out_f.write(f">{i+1} 1;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n{frag_seq}\n")
                out_f.write(f">{i+1} 2;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n{frag_seq_rc}\n")
            elif strand == '-':
                out_f.write(f">{i+1} 1;{frag_len};{tendenseq_r_fixed};{ratesseq_r_fixed}\n{frag_seq_rc}\n")
                out_f.write(f">{i+1} 2;{frag_len};{tendenseq_f_fixed};{ratesseq_f_fixed}\n{frag_seq}\n")
            i += 1
            if (i+1) % 1000000 == 0:
                t1 = time.time()
                logging.info(f"{i+1} fragments processed in {t1-t0:.2f} secs")
        if i < read_number:
            logging.warning(f"Only generated {i} fragments after {attempts} attempts. Expected {read_number}.")
    t1 = time.time()
    logging.info(f"Done processing {i} fragments in {t1-t0:.2f} secs")

### Step 7: Create reads using reseq seqToIllumina ###
def create_reads(input_fragments, reseq_model, output_reads, threads, tools):
    cmd = [tools["reseq"], "seqToIllumina", "-j", str(threads),
           "-s", reseq_model, "-i", input_fragments, "-o", output_reads, "; exit 0"]
    run_command(cmd)

### Step 8: Split interleaved FASTQ into paired gzipped FASTQ files ###
def split_reads(interleaved_fastq, output_fastq1, output_fastq2):
    """
    Split an interleaved FASTQ (each record is 4 lines) into two gzipped files:
    one for read1 and one for read2.
    """
    logging.info(f"Splitting {interleaved_fastq} into {output_fastq1} and {output_fastq2}")
    with open(interleaved_fastq, 'r') as inf, \
         gzip.open(output_fastq1, 'wt') as out1, \
         gzip.open(output_fastq2, 'wt') as out2:
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

### Step 9: Align reads using bwa mem and samtools ###
def align_reads(read1, read2, human_reference, output_bam, tools, threads=4):
    """
    Align paired-end reads with bwa mem, sort with samtools, and index the resulting BAM.
    If the bwa (or samtools) command (from the configuration) is a compound command,
    run it using bash -c.
    """
    # Process the bwa command.
    bwa_exe = tools["bwa"]
    if " " in bwa_exe:
        bwa_cmd = "bash -c '" + " ".join([bwa_exe, "mem", human_reference, read1, read2]) + "'"
        bwa_shell = True
    else:
        bwa_cmd = [bwa_exe, "mem", human_reference, read1, read2]
        bwa_shell = False

    # Process the samtools sort command.
    samtools_exe = tools["samtools"]
    if " " in samtools_exe:
        samtools_sort_cmd = "bash -c '" + " ".join([samtools_exe, "sort", "-o", output_bam, "-"]) + "'"
        samtools_shell = True
    else:
        samtools_sort_cmd = [samtools_exe, "sort", "-o", output_bam, "-"]
        samtools_shell = False

    logging.info("Aligning reads...")
    try:
        bwa_proc = subprocess.Popen(bwa_cmd, shell=bwa_shell, stdout=subprocess.PIPE, universal_newlines=True)
        sort_proc = subprocess.Popen(samtools_sort_cmd, shell=samtools_shell, stdin=bwa_proc.stdout,
                                     stdout=subprocess.PIPE, universal_newlines=True)
        bwa_proc.stdout.close()
        sort_proc.communicate()
    except Exception as e:
        logging.exception("Error during alignment")
        sys.exit(1)
    index_cmd = [tools["samtools"], "index", output_bam]
    run_command(index_cmd)

### Cleanup helper ###
def cleanup_files(file_list):
    """Remove files in the provided list if they exist."""
    for f in file_list:
        try:
            if os.path.exists(f):
                os.remove(f)
                logging.info(f"Removed intermediate file: {f}")
        except Exception as e:
            logging.warning(f"Unable to remove intermediate file {f}: {e}")

### Master orchestration function ###
def simulate_reads(config, input_fa):
    """
    Run the complete read simulation pipeline.
    
    Parameters:
      - config: a dictionary containing "tools" and "read_simulation" sections.
      - input_fa: input simulated FASTA file (e.g., muc1_simulated.fa).
    """
    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})
    
    base = os.path.splitext(input_fa)[0]
    # Intermediate files
    noNs_fa      = base + "_noNs.fasta"
    syser_fq     = base + "_noNs_syserrors.fq"
    twobit_file  = noNs_fa + ".2bit"
    subset_ref   = base + "_subset.fasta"
    psl_file     = base + "_subset.psl"
    fragments_fa = base + "_w_wessim2.fa"
    reads_fq     = base + "_w_wessim2.fq"
    # Final outputs (paired FASTQs are gzipped)
    reads_fq1    = base + "_w_wessim2_1.fq.gz"
    reads_fq2    = base + "_w_wessim2_2.fq.gz"
    output_bam   = base + "_wessim.bam"
    
    # Parameters from configuration.
    reseq_model     = rs_config.get("reseq_model")
    sample_bam      = rs_config.get("sample_bam")
    human_reference = rs_config.get("human_reference")
    read_number     = int(rs_config.get("read_number", 10000))
    fragment_size   = int(rs_config.get("fragment_size", 170))
    fragment_sd     = int(rs_config.get("fragment_sd", 35))
    min_fragment    = int(rs_config.get("min_fragment", 20))
    threads         = int(rs_config.get("threads", 24))
    bind            = int(rs_config.get("bind", 50))  # minimum fraction (%)
    
    logging.info("Starting read simulation pipeline at " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    
    # Step 1: Replace Ns.
    replace_Ns(input_fa, noNs_fa, tools)
    # Step 2: Generate systematic errors.
    generate_systematic_errors(noNs_fa, reseq_model, syser_fq, tools)
    # Step 3: Convert FASTA to 2bit.
    fa_to_twobit(noNs_fa, twobit_file, tools)
    # Step 4: Extract subset reference.
    collated_bam = extract_subset_reference(sample_bam, subset_ref, tools)
    # Step 5: Run pblat.
    run_pblat(twobit_file, subset_ref, psl_file, tools, threads=threads)
    # Step 6: Simulate fragments.
    simulate_fragments(noNs_fa, syser_fq, psl_file, read_number,
                       fragment_size, fragment_sd, min_fragment, bind, fragments_fa)
    # Check that fragments file exists and is non-empty.
    if not os.path.exists(fragments_fa) or os.path.getsize(fragments_fa) == 0:
        logging.error(f"Fragments file {fragments_fa} was not created or is empty.")
        sys.exit(1)
    # Step 7: Create reads from fragments.
    create_reads(fragments_fa, reseq_model, reads_fq, threads, tools)
    # Step 8: Split interleaved FASTQ into paired gzipped FASTQ files.
    split_reads(reads_fq, reads_fq1, reads_fq2)
    # Step 9: Align the paired reads.
    align_reads(reads_fq1, reads_fq2, human_reference, output_bam, tools, threads=threads)
    
    logging.info("Read simulation pipeline completed at " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("Final outputs:")
    logging.info("  Aligned and indexed BAM: " + output_bam)
    logging.info("  Paired FASTQ files (gzipped): " + reads_fq1 + " and " + reads_fq2)
    
    # Cleanup intermediate files.
    intermediates = [noNs_fa, syser_fq, twobit_file, subset_ref, psl_file, fragments_fa, reads_fq, collated_bam]
    cleanup_files(intermediates)

### Standalone CLI entry point ###
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

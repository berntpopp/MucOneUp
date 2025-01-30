# muc_one_up/fasta_writer.py

def write_fasta(sequences, filenames, prefix="haplotype"):
    """
    Generic function to write sequences to a FASTA file.
    `sequences`: list of strings
    `filenames`: path or file-like
    """
    with open(filenames, "w") as fh:
        for i, seq in enumerate(sequences, start=1):
            fh.write(f">{prefix}_{i}\n{seq}\n")

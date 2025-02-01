# MucOneUp

**MucOneUp** is a Python tool for simulating **MUC1 VNTR diploid references**. It builds customized references that:

1. **Generate** haplotypes containing a variable‐length VNTR region using a probability model or fixed‐lengths.  
2. **Force** a canonical terminal block (`6` or `6p` → `7 → 8 → 9`) before appending the right‐hand constant.  
3. Optionally **introduce mutations** (inserts, deletes, or replacements) in selected repeats.  
4. **Generate series of simulations** when fixed-length ranges are provided (via the `--simulate-series` flag) so that a simulation is run for each possible length (or combination of lengths for multiple haplotypes).  
5. **Run dual simulations** (normal and mutated) when a comma-separated mutation name is provided.

---

## Table of Contents

- [Installation](#installation)  
- [Quick Start](#quick-start)  
- [Usage](#usage)  
  - [Command-Line Arguments](#command-line-arguments)  
  - [Example Commands](#example-commands)  
- [Project Structure and Logic](#project-structure-and-logic)  
  - [Modules](#modules)  
  - [Config File Layout](#config-file-layout)  
- [License](#license)

---

## Installation

1. **Clone or download** this repository.  
2. **Install** in a Python 3.7+ environment:
   ```bash
   pip install .
   ```
   This will install the `muc_one_up` Python package locally.

3. (Optional) If you want to make edits and test them, install in **editable** mode:
   ```bash
   pip install -e .
   ```

Once installed, you’ll have a command-line program called **`muconeup`** available.

---

## Quick Start

1. **Create** or update a **JSON config** (see [Config File Layout](#config-file-layout)) describing your repeats, probabilities, mutations, etc.  
2. **Run** the tool by specifying your config along with desired parameters. For example:
   ```bash
   muconeup --config config.json --out-base muc1_simulated --output muc1_simulated.fa --output-structure muc1_struct.txt
   ```
3. Inspect the resulting outputs:
   - **`muc1_simulated.fa`**: multi-FASTA file of haplotype sequences.  
   - **`muc1_struct.txt`**: textual representation of each haplotype’s chain of repeats.

---

## Usage

Below are the available **command-line arguments**. Use `muconeup --help` for more details.

### Command-Line Arguments

| Argument                       | Description                                                                                                                                                                                                                                                                                                                                      |
|--------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--config <path>`              | **Required**. Path to the JSON config file containing repeats, probabilities, constants, length model, mutations, tools, and read simulation settings.                                                                                                                                                                                          |
| `--out-base <basename>`         | Base name for all output files. All outputs (simulation FASTA, VNTR structure, ORF FASTA, read simulation outputs) will be named using this base. Default is `muc1_simulated`.                                                                                                                                                           |
| `--out-dir <folder>`           | Output folder where all files will be written. Defaults to the current directory.                                                                                                                                                                                                                                                              |
| `--num-haplotypes N`           | Number of haplotypes to simulate. Typically `2` for diploid. Defaults to 2.                                                                                                                                                                                                                                                                       |
| `--fixed-lengths <vals>`       | One or more fixed lengths (or ranges) for each haplotype’s VNTR repeats. Values may be a single integer (e.g. `60`) or a range (e.g. `20-40`). When a range is provided, the default behavior is to pick one value at random from each range. Use the `--simulate-series` flag (see below) to run a simulation for every value (or combination) in the range. |
| `--simulate-series`            | (Optional) When specified and fixed-length ranges are provided, the program will generate a simulation iteration for every possible length (or combination of lengths for multiple haplotypes) instead of choosing a single random value. This flag is useful when you want to explore the entire parameter space.                             |
| `--seed <int>`                 | Random seed for reproducible simulations (affects VNTR building and mutation target selection).                                                                                                                                                                                                                                                   |
| `--mutation-name <str>`        | (Optional) Name of a mutation from the config to apply. To run dual simulations (normal and mutated), provide a comma-separated pair (e.g. `normal,dupC`). If a single value is provided, only one simulation is mutated.                                                                                                                  |
| `--mutation-targets <pairs>`   | (Optional) One or more `haplotype_index,repeat_index` pairs (1-based). E.g., `1,5 2,7`. If provided, each pair indicates which haplotype and repeat to mutate. If omitted, the mutation is applied at a random allowed repeat.                                                                                                               |
| `--output-structure`           | (Optional) If provided, output a VNTR structure file (text) listing the chain of repeats for each haplotype.                                                                                                                                                                                                                                      |
| `--output-orfs`                | (Optional) If provided, run ORF prediction and output an ORF FASTA file using the normalized naming scheme.                                                                                                                                                                                                                                     |
| `--orf-min-aa <int>`           | Minimum peptide length (in amino acids) to report from ORF prediction. Defaults to 100.                                                                                                                                                                                                                                                           |
| `--orf-aa-prefix <str>`        | (Optional) Filter resulting peptides to only those beginning with this prefix. If used without a value, defaults to `MTSSV`. If omitted, no prefix filtering is applied.                                                                                                                                                                     |
| `--simulate-reads`             | (Optional) If provided, run the read simulation pipeline on the simulated FASTA. This pipeline produces an aligned/indexed BAM and gzipped paired FASTQ files.                                                                                                                                                                              |
| `--log-level <level>`          | Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE). Default is INFO.                                                                                                                                                                                                                                                              |

### Example Commands

1. **Generate two diploid haplotypes** with random VNTR lengths (sampled from the length model defined in your config):
   ```bash
   muconeup --config config.json --out-base muc1_sim --output muc1_sim.fa --output-structure muc1_struct.txt
   ```

2. **Force a fixed length of 60 repeats** for each haplotype:
   ```bash
   muconeup --config config.json --out-base muc1_fixed --output muc1_fixed.fa --output-structure muc1_fixed.txt --num-haplotypes 2 --fixed-lengths 60
   ```

3. **Generate a single simulation using a fixed-length range** (a random value is chosen from the range for each haplotype):
   ```bash
   muconeup --config config.json --out-base muc1_random_range --fixed-lengths 20-40
   ```

4. **Generate a series of simulations** using a fixed-length range (each possible value in the range, or combination thereof, produces an output file):
   ```bash
   muconeup --config config.json --out-base muc1_series --simulate-series --fixed-lengths 20-40
   ```

5. **Apply a known mutation** (`dupC`) at a specific repeat (haplotype #1, repeat #5):
   ```bash
   muconeup --config config.json --out-base muc1_mutated --mutation-name dupC --mutation-targets 1,5
   ```

6. **Run dual simulation** (normal and mutated) with a mutation (`dupC`):
   ```bash
   muconeup --config config.json --out-base muc1_dual --mutation-name normal,dupC
   ```

7. **Apply a known mutation** (`snpA`) to a **random allowed** repeat:
   ```bash
   muconeup --config config.json --out-base muc1_random_mut --mutation-name snpA
   ```

---

## Project Structure and Logic

The **muc_one_up** Python package is organized into modules. Here is a brief summary:

```
muc_one_up/
├── cli.py           # Main CLI logic and argument parsing (now supporting series simulation and dual mutation modes)
├── config.py        # Loads and validates the JSON configuration file
├── distribution.py  # Samples the target VNTR length from a specified distribution
├── fasta_writer.py  # Helper for writing FASTA files
├── mutate.py        # Logic to apply specified mutations to haplotypes
├── probabilities.py # Provides weighted random selections for repeat transitions
├── simulate.py      # Core simulation code for building haplotypes (chains of repeats with terminal block insertion)
├── read_simulation.py  # Integrates an external read simulation pipeline to generate reads from simulated FASTA files
├── translate.py     # Translates DNA to protein and performs ORF predictions using orfipy
└── __init__.py      # Package initialization and version information
```

### High-Level Logic

1. **CLI (cli.py):**
   - Parses command-line arguments and loads the configuration.
   - If fixed-length ranges are provided, either picks a random value for each haplotype (default) or—if `--simulate-series` is specified—runs a simulation for every possible length (or combination of lengths for multiple haplotypes).
   - Simulates haplotypes via **simulate_diploid()**.
   - Optionally applies mutations using **apply_mutations()**. Dual simulation is supported when a comma-separated mutation name is provided.
   - Writes output files (FASTA, VNTR structure, ORFs) with numbered filenames.
   - Optionally runs the read simulation pipeline.

2. **Simulation (simulate.py):**
   - Constructs haplotypes by sampling repeats according to probability distributions.
   - Forces the final block of repeats (`6`/`6p` → `7` → `8` → `9`).
   - Appends left and right constant flanks from the config.

3. **Mutations (mutate.py):**
   - Applies mutations (insertion, deletion, or replacement) at specified repeats.
   - Ensures that if the current repeat symbol isn’t allowed, it is changed to an allowed one.
   - Rebuilds the haplotype sequence and marks mutated repeats with an “m” suffix.

4. **Configuration (config.py):**
   - Loads and validates the configuration JSON against a predefined schema.
   - The config file includes definitions for repeats, constants, probabilities, length model, mutations, external tool commands, and read simulation settings.

---

## Config File Layout

A simplified example:

```jsonc
{
  "repeats": {
    "1": "CACAGCATTCTTCTC...", 
    "2": "CTGAGTGGTGGAGGA...", 
    // ...
    "9": "TGAGCCTGATGCAGA..."
  },
  "constants": {
    "left": "ACGTACGTACGTACGT",
    "right": "TGCAAGCTTTGCAAGC"
  },
  "probabilities": {
    "1": { "2": 1.0 },
    "2": { "3": 1.0 },
    "3": { "4": 1.0 },
    "4": { "5": 1.0 },
    "5": { "X": 0.2, "C": 0.8 },
    // ...
    "9": { "END": 1.0 }
  },
  "length_model": {
    "distribution": "normal",
    "min_repeats": 20,
    "max_repeats": 130,
    "mean_repeats": 70,
    "median_repeats": 65
  },
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X", "C", "B", "A"], 
      "changes": [
        {
          "type": "insert",
          "start": 2,
          "end": 3,
          "sequence": "G"
        }
      ]
    }
    // Additional mutations...
  },
  "tools": {
    "reseq": "mamba run --no-capture-output -n wessim reseq",
    "faToTwoBit": "mamba run --no-capture-output -n wessim faToTwoBit",
    "samtools": "mamba run --no-capture-output -n wessim samtools",
    "pblat": "mamba run --no-capture-output -n wessim pblat",
    "bwa": "mamba run --no-capture-output -n wessim bwa"
  },
  "read_simulation": {
    "reseq_model": "reference/Hs-Nova-TruSeq.reseq",
    "sample_bam": "/path/to/sample.bam",
    "human_reference": "/path/to/GRCh38.fna.gz",
    "read_number": 10000,
    "fragment_size": 250,
    "fragment_sd": 35,
    "min_fragment": 20,
    "threads": 8
  }
}
```

---

## License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

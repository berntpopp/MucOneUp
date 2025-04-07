# MucOneUp

**MucOneUp** is a Python tool for simulating **MUC1 VNTR diploid references**. It builds customized references that:

1. **Generate** haplotypes containing a variable‐length VNTR region using a probability model or fixed‐lengths.  
2. **Force** a canonical terminal block (`6` or `6p` → `7 → 8 → 9`) before appending the right‐hand constant.  
3. Optionally **introduce mutations** (inserts, deletes, or replacements) in selected repeats with precise control over which haplotypes are affected.  
4. **Support structure files** with mutation information, allowing for reproducible generation of specific mutated VNTR structures.  
4. **Generate series of simulations** when fixed-length ranges are provided (via the `--simulate-series` flag) so that a simulation is run for each possible length (or combination of lengths for multiple haplotypes).  
5. **Run dual simulations** (normal and mutated) when a comma-separated mutation name is provided.  
6. **Detect Toxic Protein Features:** When ORF prediction is activated (via the `--output-orfs` flag), the tool scans the resulting ORF FASTA file for toxic protein sequence features by analyzing the repeat structure and amino acid composition of the variable region. A quantitative “deviation” (or similarity) score is computed relative to a wild–type model, and if the overall score exceeds a user–defined cutoff, the protein is flagged as toxic.  
7. **Generate Comprehensive Simulation Statistics:**  
   For each simulation run, detailed statistics are generated—including simulation runtime, haplotype-level metrics (repeat counts, VNTR lengths, GC content, repeat length summaries, and mutation details), as well as overall aggregated statistics. In dual simulation mode, separate reports are produced for the normal and mutated outputs.

Additionally, the **read simulation** pipeline has been integrated with a port of [w‑Wessim2](https://github.com/GeorgetteTanner/w-Wessim2). This pipeline simulates Illumina reads from the simulated FASTA by:
- Replacing Ns and generating systematic errors using external tools.
- Converting the FASTA to 2bit format.
- Extracting a subset reference from a sample BAM.
- Running pblat for alignment.
- Generating fragments and creating reads using the w‑Wessim2 port.
- Splitting interleaved FASTQ into paired FASTQ files and aligning the reads to a human reference.

---

## Table of Contents

- [Installation](#installation)
- [Installing Required Reference Files](#installing-required-reference-files)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Command-Line Arguments](#command-line-arguments)
  - [Example Commands](#example-commands)
- [Toxic Protein Detection](#toxic-protein-detection)
- [Simulation Statistics](#simulation-statistics)
- [Read Simulation Integration](#read-simulation-integration)
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

> **Optional (Conda/Mamba Environment for Read Simulation):**  
> To install all required external tools for the read simulation pipeline using conda/mamba, run:
> 
> ```bash
> mamba env create -f conda/env_wessim.yaml
> ```
> 
> After creating the environment, update the `tools` section in your configuration file to reference the executables from the newly created environment. Alternatively, you can install the tools locally.

Once installed, you’ll have a command-line program called **`muconeup`** available.

---

## Installing Required Reference Files

In order to run the read simulation pipeline, **MucOneUp** requires several external reference files (such as a human reference FASTA and a reseq model file). A helper script is provided to automate the download, extraction, and (optional) indexing of these reference files. Follow these steps to install the required references and update your configuration accordingly:

### Step 1. Run the Reference Installation Helper

A helper script is provided in the `helpers` directory. This script downloads the required reference files based on a JSON configuration (by default, `helpers/install_references_config.json`).

For example, to install the references into a directory named `./references`, run:

```bash
python helpers/install_references.py --output-dir ./references
```

This command will:
- **Download** the reference files (e.g. the GRCh38 human reference and a reseq model file).
- **Verify** the downloads via MD5 checksum.
- **Extract** any gzip-compressed files automatically.
- **Index** FASTA files using BWA if an indexing command is provided in the configuration.  
  *(If you prefer not to index the FASTA files automatically, use the `--skip-indexing` flag.)*

After the script completes, an `installed_references.json` file is generated in the output directory. This file maps each reference name to its absolute file path.

### Step 2. Update Your Configuration

Once the references are installed, update your **MucOneUp** configuration file (`config.json`) so that the fields in the `tools` and `read_simulation` sections point to the correct local paths. For example, if your reference FASTA file (e.g. GRCh38) is now located at:

```
/path/to/references/GRCh38_no_alt_analysis_set.fna.gz
```

update the `human_reference` field in the `read_simulation` section accordingly. Similarly, update the path for the reseq model file and any other reference paths you are using. For example:

```json
{
  "tools": {
    "reseq": "mamba run --no-capture-output -n env_wessim reseq",
    "faToTwoBit": "mamba run --no-capture-output -n env_wessim faToTwoBit",
    "samtools": "mamba run --no-capture-output -n env_wessim samtools",
    "pblat": "mamba run --no-capture-output -n env_wessim pblat",
    "bwa": "mamba run --no-capture-output -n env_wessim bwa"
  },
  "read_simulation": {
    "reseq_model": "/path/to/references/Hs-Nova-TruSeq.reseq",
    "sample_bam": "data/twist_v2.bam",
    "human_reference": "/path/to/references/GRCh38_no_alt_analysis_set.fna.gz",
    "read_number": 10000,
    "fragment_size": 250,
    "fragment_sd": 35,
    "min_fragment": 20,
    "threads": 8,
    "downsample_coverage": 200,
    "downsample_seed": 42,
    "reference_assembly": "hg38",
    "vntr_region_hg19": "chr1:155160500-155162000",
    "vntr_region_hg38": "chr1:155188000-155192500"
  }
}
```

Be sure to replace `/path/to/references` with the absolute path where your references were installed (as reported in the `installed_references.json` file).

### Step 3. Verify Reference Installation

Before running the full simulation or read simulation pipeline, check that the reference files are accessible and correctly referenced in your config. Review the contents of the `installed_references.json` file and ensure that each required file exists at the specified path.

---

## Quick Start

1. **Create** or update a **JSON config** (see [Config File Layout](#config-file-layout)) describing your repeats, probabilities, mutations, etc.
2. **Run** the tool by specifying your config along with desired parameters. For example:
   ```bash
   muconeup --config config.json --out-base muc1_simulated --output muc1_simulated.fa --output-structure muc1_struct.txt
   ```
3. Inspect the resulting outputs:
   - **`muc1_simulated.fa`**: Multi-FASTA file of haplotype sequences.
   - **`muc1_struct.txt`**: Textual representation of each haplotype’s chain of repeats.
   - **`simulation_stats.json`** (or with variant suffixes in dual simulation mode): JSON file(s) containing detailed simulation statistics.

---

## Usage

Below are the available **command-line arguments**. Use `muconeup --help` for more details.

### Command-Line Arguments

| Argument                       | Description                                                                                                                                                                                                                                                                                                                                      |
|--------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--config <path>`              | **Required**. Path to the JSON config file containing repeats, probabilities, constants, length model, mutations, tools, and read simulation settings.                                                                                                                                                                                          |
| `--out-base <basename>`        | Base name for all output files. All outputs (simulation FASTA, VNTR structure, ORF FASTA, read simulation outputs, and ORF toxic protein statistics) will be named using this base. Default is `muc1_simulated`.                                                                                                                             |
| `--out-dir <folder>`           | Output folder where all files will be written. Defaults to the current directory.                                                                                                                                                                                                                                                              |
| `--num-haplotypes N`           | Number of haplotypes to simulate. Typically `2` for diploid. Defaults to 2.                                                                                                                                                                                                                                                                       |
| `--fixed-lengths <vals>`       | One or more fixed lengths (or ranges) for each haplotype’s VNTR repeats. Values may be a single integer (e.g. `60`) or a range (e.g. `20-40`). When a range is provided, the default behavior is to pick one value at random from each range. Use the `--simulate-series` flag (see below) to run a simulation for every value (or combination) in the range. |
| `--simulate-series`            | (Optional) When specified and fixed-length ranges are provided, the program will generate a simulation iteration for every possible length (or combination of lengths for multiple haplotypes) instead of choosing a single random value. This flag is useful when you want to explore the entire parameter space.                             |
| `--seed <int>`                 | Random seed for reproducible simulations (affects VNTR building and mutation target selection).                                                                                                                                                                                                                                                   |
| `--mutation-name <str>`        | (Optional) Name of a mutation from the config to apply. To run dual simulations (normal and mutated), provide a comma-separated pair (e.g. `normal,dupC`). If a single value is provided, only one simulation is mutated.                                                                                                                  |
| `--mutation-targets <pairs>`   | (Optional) One or more `haplotype_index,repeat_index` pairs (1-based). E.g., `1,25 2,30`. If provided, each pair indicates which haplotype and repeat position to mutate. If omitted, the mutation is applied at a random allowed repeat. Only haplotypes specified in these targets will have mutation information in their FASTA headers.                                                                                                               |
| `--input-structure <file>`     | (Optional) Path to a structure file containing predefined VNTR repeat chains. If the structure file includes a header comment with mutation information (e.g., `# Mutation Applied: dupC (Targets: [(1, 25)])`), that information will be used to apply mutations to the specific haplotypes and positions instead of using random targets or CLI-specified targets. |
| `--output-structure`           | (Optional) If provided, output a VNTR structure file (text) listing the chain of repeats for each haplotype.                                                                                                                                                                                                                                      |
| `--output-orfs`                | (Optional) If provided, run ORF prediction and output an ORF FASTA file using the normalized naming scheme. Additionally, the resulting ORF file will be scanned for toxic protein sequence features and a JSON statistics file is generated.                                                                                                                   |
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

4. **Apply a specific mutation to targeted positions** in specific haplotypes:

   ```bash
   muconeup --config config.json --out-base muc1_mutated --mutation-name dupC --mutation-targets 1,25 2,30
   ```

   This applies the `dupC` mutation to haplotype 1 at repeat position 25 and to haplotype 2 at repeat position 30. The FASTA header for each mutated haplotype will include the mutation information.

5. **Use a structure file with mutation information** for reproducible generation of specific VNTR structures:

   ```bash
   muconeup --config config.json --out-base muc1_from_structure --input-structure structure_file.txt
   ```

   Where the structure file contains mutation information in a header comment:

   ```text
   # Mutation Applied: dupC (Targets: [(1, 25)])
   haplotype_1 1-2-3-4-5-C-X-B-X-X-X-X-X-X-X-X-V-G-B-X-X-G-A-B-Xm-X-X-X-A-A-A-B-X-D-E-C-6-7-8-9
   haplotype_2 1-2-3-4-5-C-X-A-B-X-X-X-V-G-A-A-A-B-B-X-X-X-X-X-X-X-X-X-X-X-X-X-V-V-V-V-V-V-V-V-V-G-A-B-B-X-A-A-N-R-X-X-X-X-A-B-6p-7-8-9
   ```

   Note that the mutated repeat position is marked with an "m" suffix in the structure file (e.g., `Xm`). The output FASTA will include mutation information only in the header of haplotype 1, which is the only one with a mutation according to the targets.

6. **Generate a series of simulations** using a fixed-length range (each possible value in the range, or combination thereof, produces an output file):
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

## Toxic Protein Detection

When the ORF prediction is enabled (using the `--output-orfs` flag), **MucOneUp** automatically scans the resulting ORF FASTA file for toxic protein sequence features. This new feature works as follows:

1. **Extracting the Variable Region:**  
   The ORF sequence is trimmed by removing the constant left/right flanks (if provided in the configuration), isolating the variable region (i.e. the repeat region).

2. **Detecting and Quantifying Repeats:**  
   A sliding window—whose length equals that of a consensus repeat motif (e.g. `"RCHLGPGHQAGPGLHR"`)—moves across the variable region. For each window, the similarity (using a simple Hamming distance approach) is computed, and windows exceeding a set identity threshold (e.g. 80%) are considered as valid repeats. From these, the number of repeats and the average repeat identity score are calculated.

3. **Amino Acid Composition Analysis:**  
   The tool computes the frequency of key residues (such as R, C, and H) in the variable region and compares these frequencies to a wild–type model (by default, an approximation is generated by repeating the consensus motif). A composition similarity score is calculated as:
   ```
   S_composition = 1 - (sum(|f_mut - f_wt|) / sum(f_wt))
   ```
   where a score near 1 indicates high similarity (i.e. normal) and lower scores indicate divergence (i.e. toxicity).

4. **Combining Metrics:**  
   The overall similarity (or deviation) score is computed as a weighted sum of the repeat score and the composition similarity. In this implementation, a **higher overall score indicates divergence from the wild–type** (i.e. a toxic protein), while a lower score indicates similarity to the wild–type (normal). If the overall score exceeds a user–defined toxic detection cutoff (e.g. 0.5), the ORF is flagged as toxic.

5. **Output:**  
   The detection metrics for each ORF (including repeat count, average repeat identity, repeat score, composition similarity, overall score, and the toxic flag) are saved in a JSON file (with file extension `orf_stats.txt`) alongside the ORF FASTA output.

---

## Simulation Statistics

A new feature generates a comprehensive statistics report for each simulation run. This report includes:

- **Runtime:** Total simulation runtime in seconds.
- **Haplotype-Level Metrics:** For each simulated haplotype, the report details the number of repeats, VNTR region length, GC content, individual repeat lengths (with min, max, and average), repeat type counts, and mutation details.
- **Aggregated Metrics:** Overall statistics aggregated from all haplotypes.
- **Dual Simulation Reporting:** In dual mutation mode, separate statistics reports are produced for the normal and mutated outputs.

The statistics are saved as JSON files (e.g., `muc1_simulated.002.simulation_stats.json.normal` and `muc1_simulated.002.simulation_stats.json.mut`).

---

## Read Simulation Integration

The read simulation pipeline simulates Illumina reads from the generated FASTA files. This pipeline leverages external tools (reseq, faToTwoBit, samtools, pblat, bwa) and incorporates a port of [w‑Wessim2](https://github.com/GeorgetteTanner/w-Wessim2) to:

- Replace Ns in the FASTA.
- Generate systematic errors and convert the FASTA to 2bit format.
- Extract a subset reference from a sample BAM.
- Align the 2bit file to the subset reference using pblat.
- Simulate fragments and create reads using the w‑Wessim2 port.
- Split the interleaved FASTQ into paired FASTQ files.
- Align the reads to a human reference.

---

## Structure Files with Mutation Information

MucOneUp now supports structure files that contain mutation information embedded in header comments. This allows for reproducible generation of specific mutated VNTR structures.

### Structure File Format

A structure file with mutation information has the following format:

```text
# Mutation Applied: <mutation_name> (Targets: [(<haplotype_index>, <position>), ...])
haplotype_1 <repeat_chain_with_m_suffix_for_mutated_positions>
haplotype_2 <repeat_chain>
```

Example:

```text
# Mutation Applied: dupC (Targets: [(1, 25)])
haplotype_1 1-2-3-4-5-C-X-B-X-X-X-X-X-X-X-X-V-G-B-X-X-G-A-B-Xm-X-X-X-A-A-A-B-X-D-E-C-6-7-8-9
haplotype_2 1-2-3-4-5-C-X-A-B-X-X-X-V-G-A-A-A-B-B-X-X-X-X-X-X-X-X-X-X-X-X-X-V-V-V-V-V-V-V-V-V-G-A-B-B-X-A-A-N-R-X-X-X-X-A-B-6p-7-8-9
```

When using such a structure file with the `--input-structure` argument, MucOneUp will:

1. Parse the mutation information from the header comment
2. Apply the specified mutation to the targeted haplotypes and positions
3. Generate output FASTA files with mutation information only in the headers of affected haplotypes
4. Create output structure files that preserve the mutation information

This feature is particularly useful for:

- Reproducing specific known mutations for testing
- Generating consistent test data across multiple runs
- Creating benchmark datasets for variant calling tools

## Project Structure and Logic

The **muc_one_up** Python package is organized into modules. Here is a brief summary:

```text
muc_one_up/
├── cli.py           # Main CLI logic and argument parsing (supports series simulation, dual mutation modes, toxic protein detection, simulation statistics, and read simulation)
├── config.py        # Loads and validates the JSON configuration file
├── distribution.py  # Samples the target VNTR length from a specified distribution
├── fasta_writer.py  # Helper for writing FASTA files with support for per-haplotype mutation comments
├── mutate.py        # Logic to apply specified mutations (including complex types like delete_insert) to targeted haplotypes
├── probabilities.py # Provides weighted random selections for repeat transitions
├── simulate.py      # Core simulation code for building haplotypes (chains of repeats with terminal block insertion)
├── read_simulation.py  # Integrates an external read simulation pipeline (using w‑Wessim2) to generate reads from simulated FASTA files
├── translate.py     # Translates DNA to protein and performs ORF prediction using orfipy
├── toxic_protein_detector.py   # Scans ORF FASTA outputs to detect toxic protein sequence features based on repeat structure and amino acid composition
├── simulation_statistics.py    # **New Feature**: Generates comprehensive simulation statistics for each simulation run
└── __init__.py      # Package initialization and version information
```

### High-Level Logic

1. **CLI (cli.py):**
   - Parses command-line arguments and loads the configuration.
   - Supports using predefined VNTR chains from structure files via `--input-structure`, including embedded mutation information.
   - If fixed-length ranges are provided, either picks a random value for each haplotype (default) or—if `--simulate-series` is specified—runs a simulation for every possible length (or combination of lengths for multiple haplotypes).
   - Simulates haplotypes via **simulate_diploid()** or **simulate_from_chains()** when using structure files.
   - Optionally applies mutations using **apply_mutations()** with precise control over targeted haplotypes and positions via `--mutation-targets`.
   - Dual simulation is supported when a comma-separated mutation name is provided.
   - Writes output files (FASTA, VNTR structure, ORFs) with numbered filenames and haplotype-specific mutation information in FASTA headers.
   - When ORF prediction is activated (via `--output-orfs`), the resulting ORF FASTA is further scanned for toxic protein features using **toxic_protein_detector.py**. The detection statistics are saved as a JSON file.
   - Generates a detailed simulation statistics report for each simulation iteration.
   - Optionally runs the read simulation pipeline.

2. **Simulation (simulate.py):**
   - Constructs haplotypes by sampling repeats according to probability distributions.
   - Forces the final block of repeats (`6`/`6p` → `7` → `8` → `9`).
   - Appends left and right constant flanks from the config.

3. **Mutations (mutate.py):**
   - Applies mutations (insertion, deletion, replacement, and now combined deletion–insertion via `"delete_insert"`) at specified repeats.
   - Ensures that if the current repeat symbol isn’t allowed, it is changed to an allowed one.
   - Rebuilds the haplotype sequence and marks mutated repeats with an “m” suffix.
   - Records the mutated VNTR unit sequences for separate output.

4. **Configuration (config.py):**
   - Loads and validates the configuration JSON against a predefined schema.
   - The config file includes definitions for repeats, constants, probabilities, length model, mutations, external tool commands, and read simulation settings.

---

## Config File Layout

A simplified example:

```json
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
    // Additional mutations (e.g. "delinsAT" with type "delete_insert")...
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

*Note: When using the reference installation helper (see above), update your configuration to reference the absolute paths of the downloaded files (as indicated in the `installed_references.json` file).*

---

## License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

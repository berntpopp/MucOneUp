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

Additionally, MucOneUp supports **multiple read simulation pipelines**:

**For Illumina reads**, we use a port of [w‑Wessim2](https://github.com/GeorgetteTanner/w-Wessim2). This pipeline simulates reads from the simulated FASTA by:
- Replacing Ns and generating systematic errors using external tools
- Converting the FASTA to 2bit format
- Extracting a subset reference from a sample BAM
- Running pblat for alignment
- Generating fragments and creating reads using the w‑Wessim2 port
- Splitting interleaved FASTQ into paired FASTQ files and aligning to a human reference

**For Oxford Nanopore (ONT) reads**, we integrate [NanoSim](https://github.com/bcgsc/NanoSim) to generate realistic long reads with the error profiles characteristic of nanopore sequencing. This pipeline:
- Uses pre-trained models specific to ONT technologies
- Simulates long reads with realistic error profiles
- Aligns reads to the reference using minimap2
- Generates alignment files in BAM format

---

## Table of Contents

- [Installation](#installation)
- [Installing Required Reference Files](#installing-required-reference-files)
- [Migration Guide](#migration-guide)
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

### Requirements

- Python 3.10 or higher
- pip (Python package installer)

### For Users

```bash
pip install .
```

This installs the `muconeup` command-line tool.

### For Developers

Modern Python tooling with **uv**, **ruff**, **mypy**, and automated **pre-commit hooks**.

**Setup:**
```bash
# Install uv (fast package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Initialize environment
make init
```

**Development Commands:**

| Command | Action |
|---------|--------|
| `make init` | Setup complete dev environment |
| `make test` | Run tests with coverage |
| `make lint` | Check code quality |
| `make lint-fix` | Auto-fix linting issues |
| `make format` | Format code |
| `make type-check` | Run static type checking |
| `make check` | Run all quality checks |

**Quick workflow:**
```bash
make check  # Run all checks before committing
```

> **Optional (Conda/Mamba Environments for Read Simulation):**
>
> For **Illumina read simulation** with w-Wessim2, create the environment:
> ```bash
> mamba env create -f conda/env_wessim.yaml
> ```
>
> For **ONT read simulation** with NanoSim, create the environment:
> ```bash
> mamba env create -f conda/env_nanosim.yaml
> ```
>
> After creating these environments, update the `tools` section in your configuration file to reference the executables from the newly created environments. Alternatively, you can install the tools locally.

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
    "bwa": "mamba run --no-capture-output -n env_wessim bwa",
    "nanosim": "mamba run --no-capture-output -n env_nanosim simulator.py",
    "minimap2": "mamba run --no-capture-output -n env_nanosim minimap2"
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

## Migration Guide

**Upgrading from v1.x (argparse) to v2.0 (Click)?**

MucOneUp v2.0 introduces a new Click-based CLI with improved command structure following Unix philosophy. All v1.x functionality is preserved with 100% feature parity.

**Key changes:**
- Single `muconeup` command → Four specialized commands: `simulate`, `reads`, `analyze`, `pipeline`
- Each command does ONE thing well
- Cleaner command composition

**See [docs/MIGRATION_v2.md](docs/MIGRATION_v2.md) for:**
- Complete side-by-side command comparisons
- Migration examples for all common workflows
- Breaking changes (none! - 100% backward compatible at the flag level)
- Troubleshooting tips

**Quick example:**
```bash
# v1.x (deprecated)
muconeup --config config.json --out-base sample

# v2.0 (current)
muconeup --config config.json simulate --out-base sample
```

---

## Quick Start

1. **Create** or update a **JSON config** (see [Config File Layout](#config-file-layout)) describing your repeats, probabilities, mutations, etc.

2. **Generate haplotypes:**
   ```bash
   muconeup --config config.json simulate --out-base muc1_sim --out-dir output/
   ```

3. **Inspect outputs:**
   - **`output/muc1_sim.001.simulated.fa`**: Haplotype sequences
   - **`output/muc1_sim.001.simulation_stats.json`**: Detailed statistics

4. **Optional - Add analysis:**
   ```bash
   # Predict ORFs and detect toxic proteins
   muconeup --config config.json analyze orfs output/muc1_sim.001.simulated.fa --out-base muc1_orfs

   # Generate sequence statistics
   muconeup --config config.json analyze stats output/muc1_sim.001.simulated.fa --out-base muc1_stats
   ```

5. **Optional - Simulate reads:**
   ```bash
   # Illumina short reads
   muconeup --config config.json reads illumina output/muc1_sim.001.simulated.fa --out-base muc1_reads

   # Oxford Nanopore long reads
   muconeup --config config.json reads ont output/muc1_sim.001.simulated.fa --out-base muc1_ont
   ```

---

## Usage

MucOneUp follows **Unix philosophy**: each command does one thing well. Commands can be **composed** for complete workflows.

### CLI Architecture

```bash
muconeup --config <file> [--log-level LEVEL] <command> [options]
```

**Commands:**
- **`simulate`** - Generate haplotypes ONLY (core functionality)
- **`reads`** - Simulate sequencing reads from ANY FASTA
  - `illumina` - Short read simulation
  - `ont` - Oxford Nanopore long read simulation
- **`analyze`** - Analyze ANY FASTA file
  - `orfs` - Predict ORFs and detect toxic proteins
  - `stats` - Generate sequence statistics
- **`pipeline`** - Run complete workflow (convenience orchestrator)

Use `muconeup --help` or `muconeup <command> --help` for detailed options.

### Simulate Command Options

Generate MUC1 VNTR diploid haplotypes (core functionality).

**Key Options:**
- `--out-base` - Base name for output files (default: `muc1_simulated`)
- `--out-dir` - Output directory (default: current directory)
- `--num-haplotypes` - Number of haplotypes (default: 2)
- `--seed` - Random seed for reproducibility
- `--fixed-lengths` - Fixed VNTR lengths or ranges (e.g., `60` or `20-40`)
- `--simulate-series` - Generate series across length ranges
- `--mutation-name` - Mutation to apply (or `normal,dupC` for dual mode)
- `--mutation-targets` - Specific positions to mutate (e.g., `1,25 2,30`)
- `--input-structure` - Use predefined VNTR structure file
- `--output-structure` - Write VNTR structure file
- `--random-snps` - Enable random SNP generation
- `--snp-input-file` - Apply SNPs from TSV file

**Output:**
- `{out_base}.001.simulated.fa` - Haplotype sequences
- `{out_base}.001.vntr_structure.txt` - Structure file (if `--output-structure`)
- `{out_base}.001.simulation_stats.json` - Statistics

### Reads Command Options

Simulate sequencing reads from ANY FASTA file (not just MucOneUp outputs).

**Illumina subcommand:**
```bash
muconeup --config X reads illumina INPUT_FASTA [options]
```
- `--out-base` - Base name for read files
- `--coverage` - Sequencing coverage depth
- `--threads` - Number of threads

**ONT subcommand:**
```bash
muconeup --config X reads ont INPUT_FASTA [options]
```
- `--out-base` - Base name for read files
- `--min-read-length` - Minimum read length
- `--max-read-length` - Maximum read length
- `--coverage` - Sequencing coverage depth

### Analyze Command Options

Analyze ANY FASTA file (works with external sequences).

**ORFs subcommand:**
```bash
muconeup --config X analyze orfs INPUT_FASTA [options]
```
- `--out-base` - Base name for ORF outputs
- `--orf-min-aa` - Minimum peptide length (default: 100)
- `--orf-aa-prefix` - Filter by N-terminal prefix (default: `MTSSV`)

**Stats subcommand:**
```bash
muconeup --config X analyze stats INPUT_FASTA [options]
```
- `--out-base` - Base name for statistics file

### Pipeline Command Options

Run complete workflow (convenience orchestrator).

```bash
muconeup --config X pipeline [simulate options] --with-reads [illumina|ont] --with-orfs
```

Combines `simulate` → `reads` → `analyze` in a single command.

### Example Commands

#### Basic Workflows

**1. Generate random haplotypes:**
```bash
# Generate 2 diploid haplotypes with random VNTR lengths
muconeup --config config.json simulate --out-base muc1_sim --out-dir output/
```

**2. Fixed-length haplotypes:**
```bash
# Force 60 repeats per haplotype
muconeup --config config.json simulate --out-base muc1_fixed --fixed-lengths 60

# Random pick from range 20-40
muconeup --config config.json simulate --out-base muc1_range --fixed-lengths 20-40

# Generate series (one simulation per value: 20, 21, ..., 40)
muconeup --config config.json simulate --out-base muc1_series --fixed-lengths 20-40 --simulate-series 1
```

**3. Reproducible simulations:**
```bash
# Use fixed seed for reproducibility
muconeup --config config.json simulate --out-base muc1_repro --seed 42
```

#### Mutation Workflows

**4. Apply specific mutations:**
```bash
# Mutation at specific positions (haplotype 1 repeat 25, haplotype 2 repeat 30)
muconeup --config config.json simulate --out-base muc1_mut --mutation-name dupC --mutation-targets 1,25 2,30

# Mutation at random allowed position
muconeup --config config.json simulate --out-base muc1_random_mut --mutation-name dupC
```

**5. Dual simulation (normal + mutated):**
```bash
# Generates both .normal.fa and .mut.fa outputs
muconeup --config config.json simulate --out-base muc1_dual --mutation-name normal,dupC --mutation-targets 1,25
```

**6. Use structure file:**
```bash
# Generate from predefined VNTR structure
muconeup --config config.json simulate --out-base muc1_struct --input-structure my_structure.txt
```

#### Command Composition (Step-by-Step)

**7. Complete workflow - compose commands:**
```bash
# Step 1: Generate haplotypes
muconeup --config config.json simulate --out-base muc1 --out-dir output/ --output-structure

# Step 2: Predict ORFs
muconeup --config config.json analyze orfs output/muc1.001.simulated.fa --out-base muc1_orfs

# Step 3: Generate statistics
muconeup --config config.json analyze stats output/muc1.001.simulated.fa --out-base muc1_stats

# Step 4: Simulate Illumina reads
muconeup --config config.json reads illumina output/muc1.001.simulated.fa --out-base muc1_reads --coverage 100

# Step 5: Simulate ONT reads (alternative to step 4)
muconeup --config config.json reads ont output/muc1.001.simulated.fa --out-base muc1_ont --coverage 30
```

**8. Or use pipeline for convenience:**
```bash
# All-in-one: haplotypes + ORFs + Illumina reads
muconeup --config config.json pipeline --out-base muc1_complete --out-dir output/ \
  --with-orfs --with-reads illumina

# All-in-one: haplotypes + ORFs + ONT reads
muconeup --config config.json pipeline --out-base muc1_ont_complete --out-dir output/ \
  --with-orfs --with-reads ont
```

#### SNP Integration

**9. Random SNPs:**
```bash
# Generate random SNPs (density: 1 SNP per 1000bp)
muconeup --config config.json simulate --out-base muc1_snps \
  --random-snps --random-snp-density 1.0 --random-snp-output-file snps.tsv
```

**10. Predefined SNPs:**
```bash
# Apply SNPs from file (TSV format: haplotype, position, ref_base, alt_base)
muconeup --config config.json simulate --out-base muc1_custom_snps --snp-input-file my_snps.tsv
```

**11. Dual simulation with SNPs:**
```bash
# Combine mutations and SNPs
muconeup --config config.json simulate --out-base muc1_dual_snps \
  --mutation-name normal,dupC --mutation-targets 1,25 \
  --random-snps --random-snp-density 0.5
```

#### Advanced Workflows

**12. Full pipeline with series:**
```bash
# Generate series of haplotypes (20-40 repeats) then analyze each
for i in {20..40}; do
  muconeup --config config.json simulate --out-base muc1_L${i} --fixed-lengths ${i} --seed 42
  muconeup --config config.json analyze orfs muc1_L${i}.001.simulated.fa --out-base muc1_L${i}_orfs
done
```

**13. External FASTA analysis:**
```bash
# Analyze ANY external FASTA file (not just MucOneUp outputs)
muconeup --config config.json analyze orfs /path/to/external.fa --out-base external_orfs
muconeup --config config.json reads illumina /path/to/external.fa --out-base external_reads
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

### Illumina Read Simulation

The read simulation pipeline simulates Illumina reads from the generated FASTA files. This pipeline leverages external tools (reseq, faToTwoBit, samtools, pblat, bwa) and incorporates a port of [w‑Wessim2](https://github.com/GeorgetteTanner/w-Wessim2) to:

- Replace Ns in the FASTA.
- Generate systematic errors and convert the FASTA to 2bit format.
- Extract a subset reference from a sample BAM.
- Align the 2bit file to the subset reference using pblat.
- Simulate fragments and create reads using the w‑Wessim2 port.
- Split the interleaved FASTQ into paired FASTQ files.
- Align the reads to a human reference.

---

## SNP Integration

MucOneUp supports the integration of Single Nucleotide Polymorphisms (SNPs) into simulated sequences. This feature allows for more realistic simulations by incorporating natural genetic variation. SNPs can be integrated in two ways:

1. **From a predefined file** using the `--snp-file` parameter.
2. **Generated randomly** using the `--random-snps` parameter and a specified density.

### SNP File Format

The SNP file format is tab-separated (TSV) with the following columns:

- **haplotype** (1-based): The haplotype index (1 or 2 for diploid)
- **position** (0-based): Position in the haplotype sequence
- **ref_base**: Expected reference base at the position
- **alt_base**: Alternative base to introduce

Example SNP file content:
```
haplotype	position	ref_base	alt_base
1	125	A	G
2	236	C	T
```

### Random SNP Generation

When using `--random-snps`, MucOneUp will:

1. Generate random SNPs based on the specified density (SNPs per kilobase).
2. Ensure SNPs are distributed across both haplotypes.
3. Save the generated SNPs to a file if `--random-snp-output-file` is provided.

### Dual Mutation Mode and SNPs

In dual mutation mode (using `--mutation-name normal,dupC`), SNPs are applied to both the normal and mutated sequences. For mutated sequences, the `skip_reference_check` option is automatically enabled, allowing SNPs to be applied even when mutations have altered the original reference bases.

This is particularly useful when you want to simulate scenarios where both structural mutations and SNPs are present in the sample, providing a more realistic representation of genetic diversity.

### Example Commands

```bash
# Generate random SNPs with specified density
muconeup --config config.json simulate --out-base muc1_with_snps \
  --random-snps --random-snp-density 0.5 --random-snp-output-file output/muc1_random_snps.tsv

# Apply SNPs from a predefined file
muconeup --config config.json simulate --out-base muc1_with_predefined_snps --snp-input-file my_snps.tsv

# Combine dual mutation mode with SNP integration
muconeup --config config.json simulate --out-base muc1_dual_with_snps \
  --mutation-name normal,dupC --random-snps
```

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

**Core modules:**

```text
muc_one_up/
├── cli/                    # CLI modules (config, mutations, outputs, SNPs)
├── bioinformatics/         # DNA/FASTA/SNP validation
├── read_simulator/         # Illumina and ONT pipelines
├── type_defs.py            # Type aliases and Protocol definitions
├── validation.py           # General validation functions
├── simulate.py             # Core haplotype simulation
├── mutate.py               # Mutation application
├── probabilities.py        # Weighted random selection
├── distribution.py         # Length sampling from distributions
├── fasta_writer.py         # FASTA output with mutation annotations
├── translate.py            # ORF prediction
├── toxic_protein_detector.py  # Toxic protein feature detection
└── simulation_statistics.py   # Statistics generation
```

**Key features:**
- Comprehensive type hints with mypy checking
- Runtime validation with fail-fast error messages
- 50+ test coverage with 380+ tests
- Pre-commit hooks for code quality

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

## Mutation Strict Mode

MucOneUp provides two modes for handling mutations when the target repeat is not in the mutation's `allowed_repeats` list:

### Default Mode (Auto-conversion)

By default, when a mutation is applied to a repeat that isn't listed in the `allowed_repeats` for that mutation:

1. The system automatically changes the target repeat to a randomly chosen repeat from the `allowed_repeats` list
2. A warning message is logged indicating this forced change
3. The simulation continues with the substituted repeat type

This behavior ensures simulations complete successfully even when target repeats don't match what the mutation allows.

### Strict Mode

When you need more precise control, enable strict mode by setting `"strict_mode": true` in a mutation definition:

```json
"mutations": {
  "myMutation": {
    "allowed_repeats": ["X", "C"],
    "strict_mode": true,
    "changes": [
      // mutation changes here
    ]
  }
}
```

With strict mode enabled:

1. The system validates target repeats before applying mutations
2. If a target repeat isn't in the `allowed_repeats` list, an error is raised with a detailed message
3. The simulation stops instead of automatically changing the repeat type

### Important Behavior Differences

**For explicitly specified mutation targets** (using `--mutation-targets`):

- In non-strict mode: If the target repeat isn't in `allowed_repeats`, it's automatically converted to a random allowed repeat with a warning
- In strict mode: If the target repeat isn't in `allowed_repeats`, an error is raised and the simulation stops

**For random mutation targets** (when no explicit targets provided):

- In both modes: The system only selects target positions that already have a repeat type from the `allowed_repeats` list
- This ensures that even in strict mode, randomly selected targets won't cause pipeline failures

### When to Use Strict Mode

Strict mode is particularly useful when:

- **Precision is critical**: Ensure mutations are only applied to specific repeat types
- **Debugging simulations**: Catch configuration issues early rather than having silent substitutions
- **Scientific rigor**: Prevent automatic changes that could compromise experimental design
- **Quality control**: Verify that manually specified targets meet your configuration requirements

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
      "strict_mode": false,  // Optional: When true, raises an error if target repeat isn't allowed
      "changes": [
        {
          "type": "insert",
          "start": 2,
          "end": 3,
          "sequence": "G"
        }
      ]
    },
    "delinsAT": {
      "allowed_repeats": ["C", "X"],  // Must only contain valid repeat keys from 'repeats' section
      "changes": [
        {
          "type": "delete_insert",
          "start": 2,
          "end": 4,
          "sequence": "AT"
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
  },
  "nanosim_params": {
    "training_data_path": "reference/nanosim/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1",
    "coverage": 30,
    "read_type": "ONT",
    "min_read_length": 100,
    "max_read_length": 100000,
    "threads": 8
  }
}
```

*Note: When using the reference installation helper (see above), update your configuration to reference the absolute paths of the downloaded files (as indicated in the `installed_references.json` file).*

---

## License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

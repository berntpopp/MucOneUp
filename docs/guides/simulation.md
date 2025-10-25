# VNTR Simulation Guide

Comprehensive guide to generating MUC1 VNTR diploid haplotypes with MucOneUp.

---

## Overview

The `simulate` command generates diploid haplotype sequences with customizable VNTR structures. It is the **core functionality** of MucOneUp and the starting point for all workflows.

**What it does:**

- Generates two haplotype sequences (diploid reference)
- Uses probability-based repeat transitions
- Enforces canonical terminal blocks (6/6p → 7 → 8 → 9)
- Supports fixed or random VNTR lengths
- Applies mutations (optional)
- Integrates SNPs (optional)
- Outputs FASTA, structure files, and JSON statistics

**What it does NOT do:**

- Simulate sequencing reads (use `muconeup reads` command)
- Analyze sequences (use `muconeup analyze` command)

Following Unix philosophy: **one command, one purpose**.

---

## Basic Usage

### Minimal Command

```bash
muconeup --config config.json simulate --out-base output/sample
```

**Result:** Diploid FASTA with random VNTR lengths sampled from distribution in `config.json`.

**Output:**

```
output/
└── sample.001.simulated.fa
```

---

### With Structure File

```bash
muconeup --config config.json simulate \
  --output-structure \
  --out-base output/sample
```

**Output:**

```
output/
├── sample.001.simulated.fa          # FASTA sequences
└── sample.001.vntr_structure.txt    # Repeat chain structure
```

**Structure file example:**

```
# Generated: 2025-10-20 15:30:45
# Configuration: config.json
# VNTR Lengths: haplotype_1=65 haplotype_2=58
haplotype_1 1-2-3-4-5-C-X-A-B-X-X-A-B-A-6-7-8-9
haplotype_2 1-2-3-4-5-C-X-B-X-A-X-B-A-6p-7-8-9
```

---

### With Comprehensive Statistics

```bash
muconeup --config config.json simulate \
  --output-structure \
  --output-stats \
  --out-base output/sample
```

**Output:**

```
output/
├── sample.001.simulated.fa
├── sample.001.vntr_structure.txt
└── sample.001.simulation_stats.json
```

**Statistics include:**

- Haplotype lengths, repeat counts, GC content
- Repeat chain structures
- Mutation details (if applied)
- SNP integration summary
- Runtime metrics
- Configuration snapshot

---

## VNTR Length Control

### Fixed Lengths

Generate haplotypes with exactly N repeats:

```bash
# Both haplotypes: 60 repeats
muconeup --config config.json simulate \
  --fixed-lengths 60 \
  --out-base output/fixed_60
```

**Use cases:**

- Controlled experiments
- Benchmarking at specific lengths
- Reproducible test data

---

### Random Lengths (Distribution Sampling)

Sample from distribution defined in `config.json`:

```bash
muconeup --config config.json simulate \
  --out-base output/random_length
```

**Distribution configured in config.json:**

```json
{
  "length_model": {
    "distribution_type": "normal",
    "mean_repeats": 63.3,
    "median_repeats": 70,
    "min_repeats": 42,
    "max_repeats": 85
  }
}
```

**Use cases:**

- Population-level variation
- Diverse training datasets
- Realistic heterogeneity

---

### Series Generation (Parameter Sweeps)

Generate multiple simulations across a length range:

```bash
# Generate 5 samples: 40, 50, 60, 70, 80 repeats
muconeup --config config.json simulate \
  --fixed-lengths 40-80 \
  --simulate-series 5 \
  --out-base output/series
```

**Output:**

```
output/
├── series.001.simulated.fa    # 40 repeats
├── series.002.simulated.fa    # 50 repeats
├── series.003.simulated.fa    # 60 repeats
├── series.004.simulated.fa    # 70 repeats
└── series.005.simulated.fa    # 80 repeats
```

**Progress tracking:**

```
Simulating 5 iterations  [################---]  80%  00:01:23
```

**Use cases:**

- Length-dependent benchmarking
- Coverage optimization studies
- Algorithm parameter tuning

---

## Reference Assembly Selection

Choose human genome assembly (hg19 or hg38):

```bash
# Use hg38 (default: hg19)
muconeup --config config.json simulate \
  --reference-assembly hg38 \
  --out-base output/hg38_sample
```

**Assembly differences:**

| Assembly | MUC1 Coordinates | Flanking Regions |
|----------|------------------|------------------|
| **hg19** | chr1:155,158,000-155,165,000 | Left/right constants for GRCh37 |
| **hg38** | chr1:155,185,824-155,192,916 | Left/right constants for GRCh38 |

**Important:** Ensure your reference genome matches the assembly used for simulation.

---

## Structure File Input/Output

### Generate from Predefined Structure

Use existing repeat chains instead of probability-based generation:

```bash
# Create structure file
cat > custom_structure.txt << EOF
haplotype_1 1-2-3-X-A-B-X-A-6-7-8-9
haplotype_2 1-2-X-A-X-B-A-X-6p-7-8-9
EOF

# Generate sequences from structure
muconeup --config config.json simulate \
  --input-structure custom_structure.txt \
  --out-base output/from_structure
```

**Use cases:**

- Reproduce published VNTR structures
- Test specific repeat compositions
- Validate mutation application on known structures

---

### Structure File Format

**Header:**

```
# Optional comments
# Generated: timestamp
# Configuration: config.json
```

**Body:**

```
haplotype_1 <repeat1>-<repeat2>-<repeat3>-...-<terminal_block>
haplotype_2 <repeat1>-<repeat2>-<repeat3>-...-<terminal_block>
```

**Rules:**

- Dash-separated repeat symbols
- Must end with terminal block (6 or 6p, then 7, 8, 9)
- Symbols must be defined in `config.json`
- Mutation markers (`m` suffix) added automatically when mutations applied

---

## Probability-Based Generation

### How It Works

1. **Start** with initial repeat (e.g., `1`)
2. **Sample** next repeat from probability distribution
3. **Append** to chain
4. **Repeat** until target length reached
5. **Enforce** terminal block

### Probability Matrix

Defined in `config.json`:

```json
{
  "probabilities": {
    "1": {"2": 0.3, "3": 0.2, "X": 0.5},
    "2": {"1": 0.2, "3": 0.3, "A": 0.5},
    "X": {"A": 0.4, "B": 0.4, "X": 0.2},
    "A": {"B": 0.5, "X": 0.3, "A": 0.2},
    "B": {"A": 0.4, "X": 0.4, "B": 0.2}
  }
}
```

**Example transition:**

From repeat `1`, next repeat sampled:
- `2` with 30% probability
- `3` with 20% probability
- `X` with 50% probability

### Customizing Probabilities

Update `config.json` to match real-world VNTR structures:

```bash
# Analyze existing VNTR database
muconeup --config config.json analyze vntr-stats \
  data/examples/vntr_database.tsv \
  --header \
  -o observed_probs.json

# Extract transition probabilities
jq '.transition_probabilities' observed_probs.json

# Update config.json with observed probabilities
```

---

## Output Options

### Output Directory

Specify where files are saved:

```bash
muconeup --config config.json simulate \
  --out-dir /path/to/output/ \
  --out-base sample
```

**Default:** Current working directory

---

### Output Base Name

Prefix for all output files:

```bash
muconeup --config config.json simulate \
  --out-base my_simulation
```

**Output:**

```
my_simulation.001.simulated.fa
my_simulation.001.vntr_structure.txt
my_simulation.001.simulation_stats.json
```

---

### Verbose Logging

Enable detailed logging:

```bash
muconeup --verbose --config config.json simulate \
  --out-base output/sample
```

**Alternative:**

```bash
muconeup --log-level DEBUG --config config.json simulate \
  --out-base output/sample
```

---

## Mutation Application

### Dual Simulation Mode

Generate normal and mutated pairs:

```bash
muconeup --config config.json simulate \
  --mutation-name normal,dupC \
  --mutation-targets 1,25 \
  --out-base output/dual
```

**Output:**

```
output/
├── dual.001.normal.fa                    # Normal diploid
├── dual.001.normal.vntr_structure.txt
├── dual.001.normal.simulation_stats.json
├── dual.001.mut.fa                       # Mutated diploid
├── dual.001.mut.vntr_structure.txt
└── dual.001.mut.simulation_stats.json
```

**Use case:** Benchmarking variant callers (known ground truth).

---

### Targeted Mutations

Specify exact positions:

```bash
muconeup --config config.json simulate \
  --mutation-name dupC \
  --mutation-targets 1,25 2,30 \
  --output-structure \
  --out-base output/targeted
```

**Mutation applied:**

- Haplotype 1, repeat position 25
- Haplotype 2, repeat position 30

**Structure file shows markers:**

```
# Mutation Applied: dupC (Targets: [(1, 25), (2, 30)])
haplotype_1 1-2-3-...-Xm-...-6-7-8-9
haplotype_2 1-2-3-...-Am-...-6p-7-8-9
```

The `m` suffix indicates mutated positions.

---

### Random Mutations

Let MucOneUp select positions:

```bash
muconeup --config config.json simulate \
  --mutation-name dupC \
  --random-mutation-targets 3 \
  --out-base output/random_mut
```

**Behavior:**

- Selects 3 random positions
- Respects `allowed_repeats` constraint (only mutates valid repeat types)
- Records positions in `simulation_stats.json`

---

## SNP Integration

### Random SNPs

Generate random single nucleotide polymorphisms:

```bash
muconeup --config config.json simulate \
  --random-snps \
  --random-snp-density 1.0 \
  --out-base output/with_snps
```

**Density:** SNPs per 1000 bp (1.0 = ~1 SNP per kb)

---

### File-Based SNPs

Apply SNPs from TSV file:

```bash
# Create SNP file (TSV format)
cat > snps.tsv << EOF
haplotype	position	ref	alt
1	150	A	G
1	350	C	T
2	200	G	A
EOF

# Apply SNPs
muconeup --config config.json simulate \
  --snp-input-file snps.tsv \
  --out-base output/snp_file
```

**File format:**

- **haplotype:** 1 or 2 (1-indexed)
- **position:** 0-indexed position in final sequence
- **ref:** Reference base (validated before application)
- **alt:** Alternate base

---

## Reproducibility

### Seed-Based Generation

Ensure reproducible outputs:

```bash
# Run 1
muconeup --config config.json simulate \
  --seed 42 \
  --out-base run1

# Run 2 (identical to run 1)
muconeup --config config.json simulate \
  --seed 42 \
  --out-base run2

# Verify
diff run1.001.simulated.fa run2.001.simulated.fa
# Files are identical
```

**Guarantees:**

- Same seed → identical repeat chains
- Same seed → identical random SNPs
- Platform-independent (Linux/macOS/Windows)

**Requirements:**

- Identical `config.json`
- Same MucOneUp version
- Same Python version (3.10+)

---

## Advanced Options

### Custom Configuration File

Use non-default configuration:

```bash
muconeup --config /path/to/custom_config.json simulate \
  --out-base output/sample
```

---

### Logging Levels

Control verbosity:

```bash
# Minimal output
muconeup --log-level ERROR --config config.json simulate \
  --out-base output/sample

# Detailed debug output
muconeup --log-level DEBUG --config config.json simulate \
  --out-base output/sample

# No logging
muconeup --log-level NONE --config config.json simulate \
  --out-base output/sample
```

**Levels:** DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE

---

## Common Workflows

### Benchmark Variant Caller

```bash
# Generate ground truth
muconeup --config config.json simulate \
  --mutation-name normal,dupC \
  --mutation-targets 1,25 \
  --fixed-lengths 60 \
  --output-structure \
  --out-base benchmark

# Simulate reads (separate command)
muconeup --config config.json reads illumina \
  benchmark.001.mut.fa --coverage 100
```

---

### Generate Training Dataset

```bash
# Generate 100 samples with varying lengths
for i in {1..100}; do
  muconeup --config config.json simulate \
    --seed ${i} \
    --out-base training/sample_${i}
done
```

---

### Test Mutation Detection

```bash
# Apply targeted mutation
muconeup --config config.json simulate \
  --mutation-name dupC \
  --mutation-targets 1,25 2,30 \
  --output-structure \
  --out-base mutation_test

# Check ground truth
jq '.mutations' mutation_test.001.simulation_stats.json
```

---

## Troubleshooting

### Issue: Invalid Repeat Symbol

**Error:**

```
ValueError: Invalid repeat symbol 'Z' not found in config
```

**Solution:**

Ensure all repeat symbols in structure file or mutations are defined in `config.json`:

```json
{
  "repeats": {
    "1": "AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG",
    "2": "AGTATGACCAGCAGCGTACTCTCCAGCCACAGCCCCGGTTCAGGCTCCTCCACCACTCAG",
    "Z": "SEQUENCE_FOR_Z"
  }
}
```

---

### Issue: Terminal Block Missing

**Error:**

```
ValueError: Structure missing terminal block (6/6p → 7 → 8 → 9)
```

**Solution:**

When using `--input-structure`, ensure chains end with terminal block:

```
# Correct
haplotype_1 1-2-3-X-A-B-6-7-8-9

# Incorrect (missing terminal block)
haplotype_1 1-2-3-X-A-B
```

---

### Issue: Mutation Target Invalid

**Error:**

```
MutationError: Mutation target (1, 150) exceeds haplotype length (60)
```

**Solution:**

Mutation positions must be within VNTR length:

```bash
# Haplotype has 60 repeats, so position must be ≤60
muconeup --config config.json simulate \
  --fixed-lengths 60 \
  --mutation-targets 1,25  # Valid (25 ≤ 60)
```

---

## Next Steps

**Learn More:**

- **mutations guide (coming soon)** - Apply and validate mutations
- **snps guide (coming soon)** - Advanced SNP workflows
- **[Configuration Reference](../reference/configuration.md)** - Customize repeat definitions

**Try Workflows:**

- **Workflows (coming soon)**
- **Workflows (coming soon)**

---

## Command Reference

```bash
muconeup --config CONFIG simulate [OPTIONS]

Required:
  --out-base TEXT              Output filename base

VNTR Length:
  --fixed-lengths INT|RANGE    Fixed repeat count (e.g., 60 or 40-80)
  --simulate-series INT        Number of iterations for series

Mutations:
  --mutation-name TEXT         Mutation name (or "normal,name" for dual)
  --mutation-targets TEXT      Targets as "hap,pos" (e.g., "1,25 2,30")
  --random-mutation-targets INT  Random target count

SNPs:
  --random-snps                Generate random SNPs
  --random-snp-density FLOAT   SNPs per 1000 bp (default: 1.0)
  --snp-input-file PATH        TSV file with SNPs

Structure:
  --input-structure PATH       Input structure file
  --output-structure           Output structure file

Assembly:
  --reference-assembly TEXT    hg19 or hg38 (default: hg19)

Output:
  --out-dir PATH               Output directory (default: .)
  --output-stats               Output simulation statistics JSON

Reproducibility:
  --seed INT                   Random seed for reproducibility

Logging:
  --log-level LEVEL            DEBUG|INFO|WARNING|ERROR|CRITICAL|NONE
  --verbose, -v                Enable verbose output
```

---

## See Also

- **API reference (coming soon)** - Python API for programmatic use
- **[Core Concepts](../getting-started/concepts.md)** - Understanding VNTR simulation

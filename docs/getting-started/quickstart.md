# Quick Start

Get started with MucOneUp in under 5 minutes. This tutorial walks through a complete workflow from simulation to analysis.

---

## Prerequisites

- MucOneUp installed ([Installation Guide](installation.md))
- Configuration file (`config.json` in repository root)

---

## Your First Simulation

Generate a diploid haplotype with fixed VNTR length:

```bash
# Create output directory
mkdir -p output

# Run simulation
muconeup --config config.json simulate \
  --fixed-lengths 60 \
  --out-base output/my_first_simulation
```

**What happened?**

MucOneUp generated a diploid reference with two haplotypes, each containing exactly 60 VNTR repeats.

**Output files:**

```
output/
└── my_first_simulation.001.simulated.fa
```

### View the Output

```bash
# Check file size
ls -lh output/my_first_simulation.001.simulated.fa

# View FASTA headers
grep ">" output/my_first_simulation.001.simulated.fa
# >haplotype_1 length=60_repeats assembly=hg19
# >haplotype_2 length=60_repeats assembly=hg19

# Count sequences
grep -c ">" output/my_first_simulation.001.simulated.fa
# 2  (diploid: two haplotypes)
```

---

## Add Structure Output

Generate repeat structure files to see the VNTR composition:

```bash
muconeup --config config.json simulate \
  --fixed-lengths 60 \
  --output-structure \
  --out-base output/with_structure
```

**Output files:**

```
output/
├── with_structure.001.simulated.fa
└── with_structure.001.vntr_structure.txt
```

### View the Structure

```bash
cat output/with_structure.001.vntr_structure.txt
```

**Example output:**

```
# Generated: 2025-10-20 14:45:32
# Configuration: config.json
# VNTR Length: 60 repeats
haplotype_1 1-2-3-4-5-C-X-A-B-X-X-A-B-6-7-8-9
haplotype_2 1-2-3-4-5-C-X-B-X-A-X-B-A-6p-7-8-9
```

**Understanding the structure:**

- Dash-separated symbols represent repeat units
- Terminal block (6/6p → 7 → 8 → 9) is always present
- Symbols defined in `config.json` (1, 2, 3, X, A, B, C, 6, 6p, 7, 8, 9)

---

## Apply a Mutation

Generate normal and mutated pairs to benchmark variant callers:

```bash
muconeup --config config.json simulate \
  --mutation-name normal,dupC \
  --mutation-targets 1,25 \
  --output-structure \
  --out-base output/mutation_example
```

**What happened?**

- Generated **two complete simulations**: normal and dupC mutated
- Applied dupC mutation at haplotype 1, repeat position 25
- Created structure files showing mutation markers

**Output files:**

```
output/
├── mutation_example.001.normal.fa
├── mutation_example.001.normal.vntr_structure.txt
├── mutation_example.001.normal.simulation_stats.json
├── mutation_example.001.mut.fa
├── mutation_example.001.mut.vntr_structure.txt
└── mutation_example.001.mut.simulation_stats.json
```

### View the Mutation

```bash
# Check structure file for mutation marker
cat output/mutation_example.001.mut.vntr_structure.txt
```

**Example output:**

```
# Mutation Applied: dupC (Targets: [(1, 25)])
haplotype_1 1-2-3-4-5-C-X-A-B-X-X-A-B-Xm-A-6-7-8-9
haplotype_2 1-2-3-4-5-C-X-B-X-A-X-B-A-6p-7-8-9
```

Note the `Xm` marker showing the mutated repeat position.

### Check Ground Truth

```bash
# View mutation details in JSON
jq '.mutations' output/mutation_example.001.mut.simulation_stats.json
```

**Example output:**

```json
{
  "mutation_name": "dupC",
  "targets": [[1, 25]],
  "mutated_positions": {
    "haplotype_1": [25]
  },
  "mutated_units": {
    "haplotype_1": {
      "25": "GCCCCACCCCTCCTCCCGCCGCGCCG"
    }
  }
}
```

---

## Simulate Sequencing Reads

Generate Illumina paired-end reads from your simulated haplotype:

!!! warning "Prerequisite"
    Illumina read simulation requires external tools. See [Installation: Read Simulation Setup](installation.md#read-simulation-setup).

```bash
# Simulate 30× coverage Illumina reads
muconeup --config config.json reads illumina \
  output/mutation_example.001.mut.fa \
  --coverage 30 \
  --out-base output/reads
```

**Output files:**

```
output/
├── reads_R1.fastq.gz              # Forward reads
├── reads_R2.fastq.gz              # Reverse reads
├── reads.illumina.bam             # Aligned reads
└── reads.illumina.bam.bai         # BAM index
```

### Verify Read Simulation

```bash
# Count reads
zcat output/reads_R1.fastq.gz | wc -l
# Divide by 4 to get number of reads

# Check alignment
samtools view -c output/reads.illumina.bam
# Total aligned reads

# View alignment statistics
samtools flagstat output/reads.illumina.bam
```

---

## Analyze Open Reading Frames

Predict ORFs and detect toxic protein features:

```bash
muconeup --config config.json analyze orfs \
  output/mutation_example.001.mut.fa \
  --out-base output/orfs \
  --orf-min-aa 100
```

**Output files:**

```
output/
├── orfs.pep.fa                    # Predicted peptides (FASTA)
└── orfs.orf_stats.txt             # Toxic protein statistics
```

### View ORF Statistics

```bash
cat output/orfs.orf_stats.txt
```

**Example output:**

```
Haplotype Statistics:

haplotype_1:
  Total ORFs: 12
  Toxic ORFs: 3
  Longest ORF: 245 aa
  Toxic Features:
    - ORF 3: overall_score=0.72 (TOXIC)
    - ORF 7: overall_score=0.58 (TOXIC)
    - ORF 9: overall_score=0.51 (TOXIC)

haplotype_2:
  Total ORFs: 10
  Toxic ORFs: 0
```

**Learn more:** See **[Toxic Protein Detection](../guides/toxic-protein-detection.md)** for algorithm details.

---

## Complete Workflow Example

Putting it all together: simulate, generate reads, analyze:

```bash
# 1. Generate diploid haplotype with mutation
muconeup --config config.json simulate \
  --mutation-name dupC \
  --mutation-targets 1,25 2,30 \
  --output-structure \
  --out-base output/complete_workflow

# 2. Simulate Illumina reads
muconeup --config config.json reads illumina \
  output/complete_workflow.001.simulated.fa \
  --coverage 100 \
  --out-base output/complete_reads

# 3. Predict ORFs
muconeup --config config.json analyze orfs \
  output/complete_workflow.001.simulated.fa \
  --out-base output/complete_orfs

# 4. View all outputs
ls -lh output/complete_*

# 5. Examine ground truth
jq '.' output/complete_workflow.001.simulation_stats.json
```

**What you created:**

1. Diploid reference with known dupC mutations
2. Realistic Illumina reads (100× coverage)
3. ORF predictions with toxic protein scoring
4. Ground truth data for benchmarking

**Next steps:**

- Run your variant caller on `complete_reads.illumina.bam`
- Compare variant calls to ground truth in `simulation_stats.json`
- Evaluate sensitivity and precision

---

## Batch Processing

Process multiple samples efficiently:

```bash
# Generate 10 samples with varying lengths
for i in {1..10}; do
  muconeup --config config.json simulate \
    --fixed-lengths $((40 + i * 5)) \
    --out-base output/batch_sample_${i}
done

# Simulate reads for all samples (parallel)
ls output/batch_sample_*.simulated.fa | \
  parallel muconeup --config config.json reads illumina {} \
    --coverage 30 \
    --out-base output/reads_{/.}
```

---

## Common Commands Reference

### Simulation

```bash
# Random VNTR lengths (sampled from distribution)
muconeup --config config.json simulate --out-base output/random

# Fixed VNTR lengths
muconeup --config config.json simulate \
  --fixed-lengths 60 --out-base output/fixed

# Series generation (parameter sweep)
muconeup --config config.json simulate \
  --fixed-lengths 40-80 \
  --simulate-series 5 \
  --out-base output/series

# Dual simulation (normal + mutated)
muconeup --config config.json simulate \
  --mutation-name normal,dupC \
  --out-base output/dual
```

### Read Simulation

```bash
# Illumina paired-end reads
muconeup --config config.json reads illumina \
  sample.fa --coverage 100 --out-base reads

# Oxford Nanopore long reads
muconeup --config config.json reads ont \
  sample.fa --coverage 50 --out-base reads

# PacBio HiFi reads
muconeup --config config.json reads pacbio \
  sample.fa --coverage 30 --out-base reads
```

### Analysis

```bash
# ORF prediction
muconeup --config config.json analyze orfs \
  sample.fa --out-base orfs

# Haplotype statistics
muconeup --config config.json analyze stats sample.fa

# VNTR database analysis
muconeup --config config.json analyze vntr-stats \
  database.tsv --header --structure-column vntr
```

---

## Tips and Best Practices

!!! tip "Use Structure Files"
    Always include `--output-structure` when generating test data. Structure files provide human-readable repeat chains for validation.

!!! tip "Seed for Reproducibility"
    Use `--seed` for reproducible simulations:
    ```bash
    muconeup --config config.json simulate --seed 42 --out-base reproducible
    ```

!!! tip "Start with Low Coverage"
    Test workflows with low coverage (10-30×) before running high coverage (100×+) simulations.

!!! warning "Check Configuration"
    Verify your `config.json` contains all required sections:
    - `repeats`: Repeat unit sequences
    - `constants`: Left/right flanking regions
    - `probabilities`: State transition probabilities
    - `length_model`: Distribution parameters

!!! example "Real-World Example"
    Analyze example VNTR database:
    ```bash
    muconeup --config config.json analyze vntr-stats \
      data/examples/vntr_database.tsv --header
    ```

---

## Next Steps

**Learn More:**

- **[Core Concepts](concepts.md)** - Understand VNTR simulation fundamentals
- **[Simulation Guide](../guides/simulation.md)** - Detailed simulation options
- **mutations guide (coming soon)** - Apply and validate mutations
- **read-simulation guide (coming soon)** - Platform-specific parameters

**Try Workflows:**

- **Workflows (coming soon)** - Verify pipeline sensitivity
- **Workflows (coming soon)** - Create ML datasets

---

## Getting Help

- **Documentation:** [Full Documentation](https://berntpopp.github.io/MucOneUp/)
- **Issues:** [GitHub Issue Tracker](https://github.com/berntpopp/MucOneUp/issues)
- **Examples:** `data/examples/` directory in repository

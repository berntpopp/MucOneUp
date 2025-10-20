# MucOneUp

**MUC1 VNTR simulation and analysis toolkit for genomics research**

---

## What is MucOneUp?

MucOneUp generates realistic MUC1 Variable Number Tandem Repeat (VNTR) sequences with customizable mutations and simulates sequencing reads across multiple platforms. Designed for genomics researchers studying MUC1 gene variation, it provides reproducible benchmarking, pipeline validation, and synthetic data generation for computational biology workflows.

### Why MucOneUp?

**Benchmark variant callers** with known ground truth mutations
**Test clinical pipelines** before diagnostic deployment
**Generate training data** for machine learning models
**Explore VNTR diversity** in population genetics studies
**Validate workflows** with reproducible synthetic datasets

---

## Scientific Context

The **MUC1 gene** (chr1:155,185,824-155,192,916, hg38) encodes a transmembrane glycoprotein containing a Variable Number Tandem Repeat (VNTR) region. This polymorphic region exhibits significant copy number variation (20-125 repeats) and structural complexity across human populations.

**Clinical Significance:**

- **Cancer susceptibility** - Associated with gastric, breast, and ovarian cancer risk
- **Immunological function** - Modulates innate and adaptive immune responses
- **Disease stratification** - VNTR length and mutation status correlate with disease outcomes
- **Diagnostic challenges** - Complex structure complicates accurate variant detection

Understanding MUC1 VNTR variation requires robust computational tools that generate realistic test data for algorithm development and validation.

---

## Key Features

### Realistic VNTR Simulation

Probability-based repeat transitions following biological constraints:

- Configurable repeat distributions (normal/uniform sampling)
- Canonical terminal block enforcement (6/6p → 7 → 8 → 9)
- Diploid haplotype generation with independent alleles
- Structure file support for reproducible mutations

### Flexible Mutation Engine

Insert, delete, replace, or delete-insert operations:

- Targeted or random mutation placement
- Strict validation mode (enforce allowed repeat contexts)
- Dual simulation (generate normal + mutated pairs)
- Ground truth tracking in JSON statistics

### Multi-Platform Read Simulation

Generate sequencing reads with platform-specific error profiles:

- **Illumina** - w-Wessim2 integration for paired-end reads
- **Oxford Nanopore** - NanoSim with diploid split-simulation
- **PacBio HiFi** - pbsim3 with CCS consensus accuracy
- Seed-based reproducibility across all platforms

### Comprehensive Analysis

Built-in tools for downstream analysis:

- **ORF prediction** - Identify open reading frames with orfipy
- **[Toxic protein detection](guides/toxic-protein-detection.md)** - Quantitative algorithm for ADTKD-MUC1 frameshift analysis
- **VNTR statistics** - Analyze real-world repeat databases
- **[SNaPshot validation](guides/snapshot-validation.md)** - In silico PCR → digest → extension assay simulation

### Batch Processing

Unix-style composable commands:

```bash
# Each command does one thing well
muconeup simulate --out-base sample
muconeup reads illumina sample.001.simulated.fa
muconeup analyze orfs sample.001.simulated.fa
```

Pipeline multiple files efficiently:

```bash
ls *.fa | parallel muconeup reads ont {} --coverage 50
```

---

## Quick Start

### Installation

=== "Standard Installation"

    ```bash
    # Clone repository
    git clone https://github.com/berntpopp/MucOneUp.git
    cd MucOneUp

    # Install
    make install

    # Verify installation
    muconeup --version
    ```

=== "Docker"

    ```bash
    # Pull pre-built image (includes all simulators)
    docker pull ghcr.io/berntpopp/muconeup/muconeup:latest

    # Run simulation
    docker run --rm \
      -v $(pwd)/data:/data \
      ghcr.io/berntpopp/muconeup/muconeup:latest \
      --config /app/config.json \
      simulate --out-base /data/sample
    ```

=== "Development Setup"

    ```bash
    # Modern Python tooling (uv, ruff, mypy)
    make init

    # Run tests
    make test

    # Verify code quality
    make check
    ```

### Your First Simulation

Generate diploid haplotypes with a known mutation:

```bash
# Create output directory
mkdir -p output

# Generate normal + mutated pair
muconeup --config config.json simulate \
  --mutation-name normal,dupC \
  --mutation-targets 1,25 \
  --output-structure \
  --out-base output/example

# View results
ls output/
# example.001.normal.fa
# example.001.normal.vntr_structure.txt
# example.001.normal.simulation_stats.json
# example.001.mut.fa
# example.001.mut.vntr_structure.txt
# example.001.mut.simulation_stats.json
```

**What happened?**

1. Generated two diploid references (normal and dupC mutated)
2. Applied insertion mutation at haplotype 1, repeat position 25
3. Created FASTA sequences, repeat structure files, and JSON statistics
4. Provided ground truth for benchmarking downstream tools

See [Quick Start](getting-started/quickstart.md)

---

## Research Applications

### 1. Benchmarking Variant Callers

Evaluate variant caller accuracy with known ground truth:

```bash
# Generate test data with known mutation
muconeup --config config.json simulate \
  --mutation-name normal,dupC \
  --fixed-lengths 60 \
  --out-base benchmark

# Simulate reads
muconeup --config config.json reads illumina \
  benchmark.001.mut.fa --coverage 100

# Run your variant caller
gatk HaplotypeCaller -R reference.fa \
  -I benchmark_reads.bam -O calls.vcf

# Compare to ground truth in simulation_stats.json
jq '.mutations' benchmark.001.mut.simulation_stats.json
```

**Use case:** Validate clinical pipelines before diagnostic deployment.

---

### 2. Testing Mutation Detection Pipelines

Verify your pipeline detects specific mutations:

```bash
# Apply targeted mutation
muconeup --config config.json simulate \
  --mutation-name dupC \
  --mutation-targets 1,25 2,30 \
  --output-structure \
  --out-base mutation_test

# Check if your pipeline detects the dupC mutation
# Ground truth: mutation_test.001.simulation_stats.json
```

**Dataset:** Mutation positions, affected sequences, and haplotype assignments provided in JSON.



---

### 3. Generating Synthetic Training Data

Create large-scale datasets with controlled variation:

```bash
# Generate 100 samples with varying VNTR lengths
for i in {1..100}; do
  muconeup --config config.json simulate \
    --out-base training/sample_${i} \
    --seed ${i}
done

# Simulate reads at multiple coverage levels
for cov in 30 50 100; do
  muconeup --config config.json reads illumina \
    training/*.simulated.fa --coverage ${cov}
done
```

**Use case:** Train ML models for VNTR length prediction or mutation classification.



---

### 4. Exploring Population Diversity

Analyze VNTR structures from published research:

```bash
# Analyze real-world VNTR database (44 alleles)
muconeup --config config.json analyze vntr-stats \
  data/examples/vntr_database.tsv \
  --header \
  --structure-column vntr \
  -o population_stats.json

# Extract summary statistics
jq '.mean_repeats' population_stats.json  # 63.3
jq '.median_repeats' population_stats.json  # 70
jq '.transition_probabilities' population_stats.json
```

**Dataset:** Example VNTR database with 36 unique structures (42-85 repeats).



---

### 5. Reproducible Research

Ensure reproducibility with seed-based generation:

```bash
# Same seed → identical output (deterministic)
muconeup --config config.json simulate \
  --seed 42 \
  --out-base reproducible

muconeup --config config.json reads illumina \
  reproducible.001.simulated.fa \
  --seed 42 \
  --out-base reads

# Share config + seeds → fully reproducible datasets
```

**Publication-ready:** Enables method comparisons and collaborative research.



---

## Example Output

MucOneUp generates comprehensive outputs for each simulation:

```
output/
├── sample.001.simulated.fa          # FASTA sequences (diploid haplotypes)
├── sample.001.vntr_structure.txt    # Repeat chain (1-2-3-X-A-B-6-7-8-9)
├── sample.001.simulation_stats.json # Metrics and ground truth
├── sample.pep.fa                     # Predicted ORF peptides
├── sample.orf_stats.txt              # Toxic protein scores
└── sample.illumina.bam               # Simulated reads (aligned)
```

**Statistics include:**

- Haplotype lengths, repeat counts, GC content
- Mutation positions and affected sequences
- SNP integration summary with validation results
- Runtime metrics and configuration snapshot

---

## Documentation

<div class="grid cards" markdown>

-  **[Getting Started](getting-started/installation.md)**
  Installation, quick start tutorial, and core concepts

-  **[User Guides](guides/simulation.md)**
  Detailed documentation for simulation, mutations, SNPs, and analysis

-  **[Workflows](workflows/)**
  Real-world research applications with complete examples

-  **[CLI Reference](reference/cli.md)**
  Complete command-line interface documentation

-  **API Reference (coming soon)**
  Python API documentation auto-generated from source

  Analysis of real-world VNTR databases

  Development setup, architecture, testing, and code style

-  **[About](about/citation.md)**
  Citation guide, license, and changelog

</div>

---

## Citation

If you use MucOneUp in your research, please cite:

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  url = {https://github.com/berntpopp/MucOneUp},
  note = {Software version available at https://github.com/berntpopp/MucOneUp/releases}
}
```

A manuscript describing MucOneUp is in preparation.

See [Citation Guide](about/citation.md)

---

## Community

**GitHub Repository:** [berntpopp/MucOneUp](https://github.com/berntpopp/MucOneUp)
**Docker Images:** [GitHub Container Registry](https://github.com/berntpopp/MucOneUp/pkgs/container/muconeup%2Fmuconeup)
**Issue Tracker:** [Report bugs or request features](https://github.com/berntpopp/MucOneUp/issues)
**Discussions:** Questions and community support

[![GitHub Stars](https://img.shields.io/github/stars/berntpopp/MucOneUp?style=social)](https://github.com/berntpopp/MucOneUp)
[![Tests](https://github.com/berntpopp/MucOneUp/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneUp/actions)
[![Docker](https://img.shields.io/badge/docker-ghcr.io-blue.svg)](https://github.com/berntpopp/MucOneUp/pkgs/container/muconeup%2Fmuconeup)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

**Development Status:** Active | **License:** MIT | **Maintained by:** [Bernt Popp](https://github.com/berntpopp)

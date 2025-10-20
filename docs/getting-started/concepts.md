# Core Concepts

Understanding the fundamental concepts behind MucOneUp helps you design better experiments and interpret results accurately.

---

## MUC1 Gene and VNTR Region

### Biological Background

**MUC1** (Mucin 1) is a transmembrane glycoprotein encoded on chromosome 1 (chr1:155,185,824-155,192,916, hg38). It plays critical roles in:

- **Epithelial protection** - Physical barrier on cell surfaces
- **Cell signaling** - Regulates cell growth and differentiation
- **Immune modulation** - Innate and adaptive immune responses

### Variable Number Tandem Repeat (VNTR)

The MUC1 gene contains a **polymorphic VNTR region** consisting of tandem repeats of a 60-bp consensus sequence. This region exhibits:

**Structural Characteristics:**

- **Repeat unit:** 60 bp encoding 20 amino acids (PAHGVTSAPDTRPAPGSTAPPA)
- **Copy number variation:** 20-125 repeats across human populations
- **Length polymorphism:** Different alleles vary in repeat count
- **Sequence variation:** Individual repeats show nucleotide polymorphisms

**Clinical Significance:**

- **Cancer association** - Altered VNTR length correlates with cancer risk (gastric, breast, ovarian)
- **Immune function** - VNTR structure affects antigen presentation
- **Diagnostic challenges** - Complex structure complicates accurate sequencing and variant calling

### Why Simulate MUC1 VNTR?

Real MUC1 VNTR sequencing presents challenges:

- **Alignment ambiguity** - Repetitive structure causes mapping errors
- **Variant calling difficulty** - Standard pipelines struggle with tandem repeats
- **Lack of ground truth** - Real data doesn't provide known mutation positions

**MucOneUp solves this by:**

- Generating realistic VNTR sequences with biological constraints
- Providing ground truth for mutations and structural variants
- Enabling controlled benchmarking of bioinformatics tools

---

## Diploid Haplotypes

### What are Diploid Haplotypes?

Humans are **diploid organisms** - we inherit two copies of each chromosome (one maternal, one paternal). Therefore, each individual has **two MUC1 VNTR alleles**.

**MucOneUp generates diploid references:**

```
Individual Genotype:
├── Haplotype 1 (Maternal): 1-2-3-X-A-B-6-7-8-9  (60 repeats)
└── Haplotype 2 (Paternal): 1-2-X-A-X-B-6p-7-8-9 (58 repeats)
```

### Why Diploid Matters

**Variant Calling:**

- Variant callers expect diploid genomes (heterozygous vs homozygous calls)
- Simulating only one haplotype produces unrealistic data
- Diploid simulation enables proper benchmarking

**Read Simulation:**

- Reads sample from both haplotypes (allelic balance)
- Coverage distributes across maternal and paternal alleles
- Mimics real sequencing experiments

**Mutation Testing:**

- Test heterozygous mutations (one haplotype affected)
- Test homozygous mutations (both haplotypes affected)
- Evaluate phasing accuracy

---

## VNTR Structure and Repeat Units

### Repeat Symbols

MucOneUp represents VNTR structure using symbolic notation defined in `config.json`:

**Core Repeats:**

| Symbol | Description | Example Sequence (first 30 bp) |
|--------|-------------|------------------|
| `1` | Canonical repeat variant 1 | AAGGAGACTTCGGCTACCCAGAGAAG... |
| `2` | Canonical repeat variant 2 | AGTATGACCAGCAGCGTACTCTCCAG... |
| `3` | Canonical repeat variant 3 | GGACAGGATGTCACTCTGGCCCCGGC... |
| `X` | Variable repeat | GCCCACGGTGTCACCTCGGCCCCGGA... |
| `A` | Polymorphism type A | GCCCACGGTGTCACCTCGGCCCCGGA... |
| `B` | Polymorphism type B | GCCCACGGTGTCACCTCGGCCCCGGA... |
| `C` | Polymorphism type C | GCCCACGGTGTCACCTCGGCCCCGGA... |

**Terminal Block (Always Present):**

| Symbol | Description |
|--------|-------------|
| `6` | Pre-terminal repeat variant 1 |
| `6p` | Pre-terminal repeat variant 2 (polymorphic) |
| `7` | Terminal repeat 1 |
| `8` | Terminal repeat 2 |
| `9` | Terminal repeat 3 (final) |

### Canonical Terminal Block

MucOneUp **enforces** a canonical terminal block matching biological MUC1 structure:

```
Variable Region → 6 or 6p → 7 → 8 → 9 → Constant Region
```

**Why this matters:**

- Biological MUC1 always ends with this conserved block
- Simulating without it produces unrealistic sequences
- Terminal block prevents premature chain termination

**Example:**

```
# Realistic (enforced by MucOneUp)
haplotype_1: 1-2-3-X-A-B-X-6-7-8-9

# Unrealistic (not allowed)
haplotype_1: 1-2-3-X-A-B-X-9  # Missing 6-7-8
```

---

## Probability-Based Repeat Selection

### How Repeat Chains are Generated

MucOneUp uses **state transition probabilities** to generate realistic repeat chains:

**Process:**

1. **Start** with initial repeat (e.g., `1`)
2. **Sample** next repeat from probability distribution
3. **Append** to chain
4. **Repeat** until target length reached
5. **Enforce** terminal block (6/6p → 7 → 8 → 9)

### Probability Matrix

Defined in `config.json`:

```json
{
  "probabilities": {
    "1": {"2": 0.3, "3": 0.2, "X": 0.5},
    "2": {"1": 0.2, "3": 0.3, "A": 0.5},
    "X": {"A": 0.4, "B": 0.4, "X": 0.2}
  }
}
```

**Interpretation:**

- From repeat `1`, transition to:
  - `2` with 30% probability
  - `3` with 20% probability
  - `X` with 50% probability

### Why Probability-Based?

**Biological Realism:**

- Real VNTR sequences show non-random repeat patterns
- Certain transitions occur more frequently than others
- Probability model captures biological constraints

**Reproducibility:**

- With same seed, generates identical sequences
- Enables controlled experiments

**Customization:**

- Update probabilities to match specific populations
- Simulate rare or common haplotype structures

---

## VNTR Length Sampling

### Length Distributions

MucOneUp supports two distribution models:

**Normal Distribution:**

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

Samples repeat counts from normal distribution, clipped to [min, max] range.

**Uniform Distribution:**

```json
{
  "length_model": {
    "distribution_type": "uniform",
    "min_repeats": 40,
    "max_repeats": 100
  }
}
```

Samples uniformly across [min, max] range.

### Fixed vs Random Lengths

**Fixed Lengths:**

```bash
# Both haplotypes get exactly 60 repeats
muconeup --config config.json simulate --fixed-lengths 60
```

**Random Lengths:**

```bash
# Sample from distribution in config.json
muconeup --config config.json simulate
```

**When to use fixed:**

- Controlled experiments
- Benchmarking across specific lengths
- Reproducible test cases

**When to use random:**

- Population-level variation
- Training data diversity
- Realistic heterogeneity

---

## Mutations

### Mutation Types

MucOneUp supports four mutation operations:

**1. Insert**

Add sequence at specific position:

```
Before: 1-2-3-X-A-B-6-7-8-9
After:  1-2-3-X-INSERTED-A-B-6-7-8-9
```

**2. Delete**

Remove repeat at specific position:

```
Before: 1-2-3-X-A-B-6-7-8-9
After:  1-2-3-A-B-6-7-8-9  (X deleted)
```

**3. Replace**

Substitute repeat at specific position:

```
Before: 1-2-3-X-A-B-6-7-8-9
After:  1-2-3-REPLACED-A-B-6-7-8-9
```

**4. Delete-Insert**

Combine delete and insert operations:

```
Before: 1-2-3-X-A-B-6-7-8-9
After:  1-2-3-NEW_SEQ-B-6-7-8-9  (X deleted, NEW_SEQ inserted)
```

### Mutation Targeting

**Targeted Mutations:**

Specify exact haplotype and position:

```bash
muconeup --config config.json simulate \
  --mutation-name dupC \
  --mutation-targets 1,25 2,30
```

Applies dupC mutation to:
- Haplotype 1, position 25
- Haplotype 2, position 30

**Random Mutations:**

Let MucOneUp select positions:

```bash
muconeup --config config.json simulate \
  --mutation-name dupC \
  --random-mutation-targets 2
```

Randomly selects 2 positions respecting `allowed_repeats` constraints.

### Mutation Validation

**Allowed Repeats:**

Mutations define which repeat contexts are valid:

```json
{
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X", "A", "B"],
      "changes": [...]
    }
  }
}
```

dupC mutation can only be applied to X, A, or B repeats.

**Strict Mode:**

```json
{
  "mutations": {
    "dupC": {
      "strict_mode": true,
      "allowed_repeats": ["X"]
    }
  }
}
```

- **strict_mode: true** - Error if target repeat not in allowed_repeats
- **strict_mode: false** - Auto-convert to nearest allowed repeat (with warning)

---

## Ground Truth and Statistics

### Simulation Statistics

Every simulation generates a JSON file with comprehensive metadata:

```json
{
  "timestamp": "2025-10-20T14:45:32",
  "configuration": {
    "config_file": "config.json",
    "reference_assembly": "hg19"
  },
  "haplotypes": {
    "haplotype_1": {
      "sequence_length": 12450,
      "repeat_count": 60,
      "gc_content": 0.58,
      "repeat_chain": "1-2-3-X-A-B-..."
    },
    "haplotype_2": { ... }
  },
  "mutations": {
    "mutation_name": "dupC",
    "targets": [[1, 25]],
    "mutated_positions": { ... },
    "mutated_units": { ... }
  },
  "runtime_seconds": 2.34
}
```

### Using Ground Truth

**Benchmarking Variant Callers:**

1. Extract mutation positions from `simulation_stats.json`
2. Run variant caller on simulated reads
3. Compare caller VCF to ground truth
4. Calculate sensitivity (true positive rate)

**Example:**

```bash
# Extract ground truth
jq '.mutations.targets' sample.simulation_stats.json
# [[1, 25]]

# Check if variant caller detected position 25 on haplotype 1
bcftools view calls.vcf | grep -A5 "POS.*25"
```

---

## SNP Integration

### What are SNPs?

Single Nucleotide Polymorphisms (SNPs) are single-base variations in DNA sequence. MucOneUp integrates SNPs into haplotypes for increased realism.

### SNP Application Methods

**Random SNPs:**

```bash
muconeup --config config.json simulate \
  --random-snps \
  --random-snp-density 1.0 \
  --out-base snp_sim
```

Generates ~1 SNP per 1000 bp (density = 1.0).

**File-Based SNPs:**

```bash
muconeup --config config.json simulate \
  --snp-input-file variants.tsv \
  --out-base snp_sim
```

**SNP File Format (TSV):**

```
haplotype	position	ref	alt
1	150	A	G
1	350	C	T
2	200	G	A
```

- **haplotype:** 1 or 2 (1-indexed)
- **position:** 0-indexed position in final sequence
- **ref:** Reference base (validated before application)
- **alt:** Alternate base

### SNP Validation

MucOneUp validates reference bases before applying SNPs:

```
Position 150: Expected ref=A, Found=A → Applied (A→G)
Position 350: Expected ref=C, Found=G → Skipped (reference mismatch)
```

**Statistics report:**

```json
{
  "snp_integration": {
    "attempted": 10,
    "successful": 8,
    "failed": 2,
    "failure_reasons": {
      "reference_mismatch": 2
    }
  }
}
```

---

## Read Simulation

### Why Simulate Reads?

Real sequencing data contains:

- **Platform-specific errors** - Illumina substitution bias, ONT homopolymer errors
- **Coverage variation** - Uneven depth across genome
- **Quality scores** - Per-base confidence values
- **Fragment lengths** - Insert size distributions (paired-end)

**Simulating reads enables:**

- Testing alignment algorithms
- Benchmarking variant callers with realistic error profiles
- Evaluating coverage requirements
- Training sequencing-aware ML models

### Platform Differences

**Illumina (Short Reads):**

- **Read length:** 100-300 bp (paired-end)
- **Error rate:** 0.1-1%
- **Error type:** Substitutions (A↔C bias in chemistry)
- **Coverage:** High (50-150×)
- **Use case:** SNV detection, short indels

**Oxford Nanopore (Long Reads):**

- **Read length:** 1-100 kb
- **Error rate:** 5-15% (raw), 1-5% (consensus)
- **Error type:** Insertions/deletions (homopolymer errors)
- **Coverage:** Moderate (30-50×)
- **Use case:** Structural variants, phasing, repetitive regions

**PacBio HiFi (Long Accurate Reads):**

- **Read length:** 10-25 kb
- **Error rate:** 0.1-1% (CCS consensus)
- **Error type:** Random (no systematic bias)
- **Coverage:** Moderate (30-50×)
- **Use case:** Structural variants, phasing, high accuracy

### Diploid Split-Simulation (ONT/PacBio)

**Challenge:** Long-read simulators sample reads proportional to sequence length. In diploid references, longer haplotypes receive disproportionately more reads (allelic bias).

**MucOneUp solution:**

1. **Detect** diploid reference (exactly 2 sequences)
2. **Split** into separate haplotype files
3. **Simulate** each haplotype independently (equal coverage)
4. **Merge** reads from both haplotypes
5. **Align** merged reads to diploid reference

**Result:** Balanced allelic coverage.

---

## Reproducibility

### Seed-Based Determinism

Use `--seed` for reproducible outputs:

```bash
# Run 1
muconeup --config config.json simulate --seed 42 --out-base run1
muconeup --config config.json reads illumina run1.001.simulated.fa --seed 42

# Run 2 (identical to run 1)
muconeup --config config.json simulate --seed 42 --out-base run2
muconeup --config config.json reads illumina run2.001.simulated.fa --seed 42
```

**Guarantees:**

- Same seed → identical VNTR structures
- Same seed → identical read sampling
- Platform-independent (Linux/macOS/Windows)

**Requirements:**

- Identical `config.json`
- Same tool versions (MucOneUp, reseq, NanoSim, etc.)
- Same Python version (3.10+)

### Sharing Reproducible Datasets

**For publications:**

1. Share `config.json`
2. Document seeds used
3. Specify tool versions
4. Provide structure files (human-readable validation)

**Example:**

```
Methods:
Synthetic MUC1 VNTR sequences generated with MucOneUp v0.19.0
(seed=42, config.json provided in supplementary materials).
Illumina reads simulated at 100× coverage (reseq v1.1, seed=42).
```

---

## Next Steps

**Understand the Tools:**

- **[Simulation Guide](../guides/simulation.md)** - Detailed VNTR generation options
- **mutations guide (coming soon)** - Apply and validate mutations
- **read-simulation guide (coming soon)** - Platform-specific parameters

**Try Workflows:**

- **Workflows (coming soon)** - Seed-based workflows

**Reference:**

- **[Configuration Guide](../reference/configuration.md)** - Customize repeat definitions and probabilities
- **API reference (coming soon)** - Python API for custom workflows

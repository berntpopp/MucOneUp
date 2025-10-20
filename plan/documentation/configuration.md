# Configuration Reference

Complete reference for `config.json` structure and customization.

---

## Overview

MucOneUp requires a JSON configuration file specifying:

- **Repeat definitions** - DNA sequences for each repeat symbol
- **Flanking constants** - Left/right regions surrounding VNTR
- **Probability transitions** - State transitions for repeat selection
- **Length model** - Distribution parameters for VNTR lengths
- **Mutation definitions** - Named mutations with operations
- **Tool paths** - External tools for read simulation
- **Simulation parameters** - Platform-specific read simulation settings

**Location:** Repository root contains `config.json` (example configuration).

---

## File Structure

```json
{
  "repeats": { ... },
  "constants": { ... },
  "probabilities": { ... },
  "length_model": { ... },
  "mutations": { ... },
  "tools": { ... },
  "read_simulation": { ... },
  "nanosim_params": { ... },
  "pacbio_params": { ... }
}
```

---

## Repeats Section

### Purpose

Defines DNA sequences for each repeat symbol used in VNTR chains.

### Structure

```json
{
  "repeats": {
    "1": "GCCCCACCCCTCCTCCCGCCGCGCCG",
    "2": "GCTCCACCCCTTCTCCCACCGCGCCG",
    "3": "GCCCCACCCCTTCTCCCACCGCGCCG",
    "X": "GCCCCACCCCTCCTCCCGCCGCGCCG",
    "A": "GCTCCACCCCTTCTCCCACCGCGCCG",
    "B": "GCCCCACCCCTTCTCCCACCGCGCCG",
    "C": "GCCCCACCCCTCTTCCCGCCGCGCCG",
    "6": "GCCCCACCCCTCCTCCCGCCGCGCCG",
    "6p": "GCTCCACCCCTTCTCCCACCGCGCCG",
    "7": "GCCCCACCCCTTCTCCCACCGCGCCG",
    "8": "GCCCCACCCCTCCTCCCGCCGCGCCG",
    "9": "GCTCCACCCCTTCTCCCACCGCGCCG"
  }
}
```

### Rules

- **Key:** Repeat symbol (string, typically 1-2 characters)
- **Value:** DNA sequence (uppercase ACGT)
- **Terminal block required:** Must define 6 or 6p, 7, 8, 9
- **Validation:** Sequences validated on load (must be valid DNA)

### Example: Custom Repeat

```json
{
  "repeats": {
    "CUSTOM": "ATCGATCGATCGATCGATCGATCGAT"
  },
  "probabilities": {
    "1": {"CUSTOM": 0.2, "2": 0.8}
  }
}
```

---

## Constants Section

### Purpose

Defines left and right flanking regions surrounding the VNTR, assembly-specific.

### Structure (Nested Format)

```json
{
  "constants": {
    "hg19": {
      "left": "AGCAGGCAGTGCGGGGCCGCTGCTGCTG...",
      "right": "TGCTGCTGCTGCGGGGCCGCTGCTGCTG...",
      "vntr_start": 1000,
      "vntr_end": 8000
    },
    "hg38": {
      "left": "AGCAGGCAGTGCGGGGCCGCTGCTGCTG...",
      "right": "TGCTGCTGCTGCGGGGCCGCTGCTGCTG...",
      "vntr_start": 1500,
      "vntr_end": 8500
    }
  }
}
```

### Structure (Flat Format, Auto-Converted)

```json
{
  "constants": {
    "left": "AGCAGGCAGTGCGGGGCCGCTGCTGCTG...",
    "right": "TGCTGCTGCTGCGGGGCCGCTGCTGCTG...",
    "vntr_start": 1000,
    "vntr_end": 8000
  }
}
```

**Note:** Flat format assumed to be hg19 and auto-converted to nested format.

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `left` | string | DNA sequence upstream of VNTR |
| `right` | string | DNA sequence downstream of VNTR |
| `vntr_start` | integer | VNTR start position (0-indexed) |
| `vntr_end` | integer | VNTR end position (0-indexed) |

### Assembly-Specific Constants

MUC1 genomic coordinates differ between assemblies:

| Assembly | Chromosome | Start | End |
|----------|------------|-------|-----|
| **hg19** | chr1 | 155,158,000 | 155,165,000 |
| **hg38** | chr1 | 155,185,824 | 155,192,916 |

Ensure `left` and `right` constants match your chosen assembly.

---

## Probabilities Section

### Purpose

State transition probabilities for repeat selection during chain generation.

### Structure

```json
{
  "probabilities": {
    "1": {"2": 0.3, "3": 0.2, "X": 0.5},
    "2": {"1": 0.2, "3": 0.3, "A": 0.5},
    "3": {"1": 0.1, "2": 0.4, "X": 0.5},
    "X": {"A": 0.4, "B": 0.4, "X": 0.2},
    "A": {"B": 0.5, "X": 0.3, "A": 0.2},
    "B": {"A": 0.4, "X": 0.4, "B": 0.2},
    "C": {"X": 0.6, "A": 0.4}
  }
}
```

### Rules

- **Outer key:** Source repeat symbol
- **Inner key:** Destination repeat symbol
- **Value:** Transition probability (0.0 to 1.0)
- **Sum:** Inner values should sum to 1.0 (normalized automatically if not)

### Example

From repeat `1`, transitions:

- To `2` with 30% probability
- To `3` with 20% probability
- To `X` with 50% probability

### Customizing Probabilities

**From Real Data:**

```bash
# Analyze VNTR database
muconeup --config config.json analyze vntr-stats \
  data/examples/vntr_database.tsv \
  --header \
  -o observed_probs.json

# Extract transitions
jq '.transition_probabilities' observed_probs.json

# Update config.json
```

**Manual Specification:**

```json
{
  "probabilities": {
    "X": {"X": 0.8, "A": 0.2}  # X repeats frequently
  }
}
```

---

## Length Model Section

### Purpose

Distribution parameters for sampling VNTR repeat counts.

### Normal Distribution

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

**Behavior:**

- Sample from normal distribution (μ=63.3, σ derived from median)
- Clip to [42, 85] range

### Uniform Distribution

```json
{
  "length_model": {
    "distribution_type": "uniform",
    "min_repeats": 40,
    "max_repeats": 100
  }
}
```

**Behavior:**

- Sample uniformly from [40, 100]

### Fields

| Field | Type | Description | Required |
|-------|------|-------------|----------|
| `distribution_type` | string | "normal" or "uniform" | Yes |
| `mean_repeats` | float | Mean for normal distribution | If normal |
| `median_repeats` | float | Median for normal distribution | If normal |
| `min_repeats` | integer | Minimum repeat count | Yes |
| `max_repeats` | integer | Maximum repeat count | Yes |

---

## Mutations Section

### Purpose

Define named mutations with operations (insert/delete/replace/delete_insert).

### Structure

```json
{
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X", "A", "B"],
      "strict_mode": false,
      "changes": [
        {
          "operation": "insert",
          "sequence": "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA",
          "position_offset": 0
        }
      ]
    },
    "deletion_example": {
      "allowed_repeats": ["X"],
      "strict_mode": true,
      "changes": [
        {
          "operation": "delete"
        }
      ]
    }
  }
}
```

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `allowed_repeats` | array | Valid repeat symbols for this mutation |
| `strict_mode` | boolean | Enforce allowed_repeats (error if violated) |
| `changes` | array | List of mutation operations |

### Mutation Operations

**Insert:**

```json
{
  "operation": "insert",
  "sequence": "ATCGATCGATCG",
  "position_offset": 0
}
```

Inserts sequence at target position.

**Delete:**

```json
{
  "operation": "delete"
}
```

Removes repeat at target position.

**Replace:**

```json
{
  "operation": "replace",
  "sequence": "ATCGATCGATCG"
}
```

Substitutes repeat at target position with new sequence.

**Delete-Insert:**

```json
{
  "operation": "delete_insert",
  "sequence": "ATCGATCGATCG",
  "position_offset": 0
}
```

Deletes repeat, then inserts sequence.

### Strict Mode

**strict_mode: false (Permissive):**

- If target repeat not in `allowed_repeats`, auto-convert to nearest allowed repeat
- Emit warning
- Simulation continues

**strict_mode: true (Strict):**

- If target repeat not in `allowed_repeats`, raise error
- Simulation fails

**Example:**

```json
{
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X"],
      "strict_mode": true,
      "changes": [...]
    }
  }
}
```

```bash
# This will fail if position 25 is not an "X" repeat
muconeup --config config.json simulate \
  --mutation-name dupC \
  --mutation-targets 1,25
```

---

## Tools Section

### Purpose

Paths to external tools for read simulation.

### Structure

```json
{
  "tools": {
    "reseq": "/path/to/reseq",
    "bwa": "/usr/bin/bwa",
    "samtools": "/usr/bin/samtools",
    "faToTwoBit": "/usr/bin/faToTwoBit",
    "pblat": "/usr/bin/pblat",
    "minimap2": "/usr/bin/minimap2",
    "pbsim3": "/usr/bin/pbsim3",
    "ccs": "/usr/bin/ccs"
  }
}
```

### Auto-Detection

If paths not specified, MucOneUp searches system PATH:

```json
{
  "tools": {}  // Auto-detect all tools
}
```

### Conda Environments

When using conda environments, specify full paths:

```bash
# Find tool path in conda env
conda activate wessim
which reseq
# /home/user/miniconda3/envs/wessim/bin/reseq

# Update config.json
{
  "tools": {
    "reseq": "/home/user/miniconda3/envs/wessim/bin/reseq"
  }
}
```

---

## Read Simulation Section

### Purpose

Parameters for Illumina read simulation (w-Wessim2 pipeline).

### Structure

```json
{
  "read_simulation": {
    "simulator": "illumina",
    "read_length": 150,
    "fragment_size": 350,
    "fragment_sd": 50,
    "coverage": 100,
    "threads": 4,
    "reference_genome": "/path/to/hg38.fa",
    "error_model": "reseq_illumina",
    "seed": null
  }
}
```

### Fields

| Field | Type | Description | Default |
|-------|------|-------------|---------|
| `simulator` | string | "illumina", "ont", or "pacbio" | "illumina" |
| `read_length` | integer | Read length (bp) | 150 |
| `fragment_size` | integer | Mean insert size (bp) | 350 |
| `fragment_sd` | integer | Insert size std dev (bp) | 50 |
| `coverage` | integer | Target coverage depth | 100 |
| `threads` | integer | Parallel threads | 4 |
| `reference_genome` | string | Path to reference FASTA | Required |
| `error_model` | string | Error model name | "reseq_illumina" |
| `seed` | integer | Random seed (null = random) | null |

---

## NanoSim Parameters Section

### Purpose

Parameters for Oxford Nanopore read simulation.

### Structure

```json
{
  "nanosim_params": {
    "training_data_path": "/path/to/nanosim/training",
    "coverage": 50,
    "min_read_length": 1000,
    "max_read_length": 10000,
    "correction_factor": 0.325,
    "enable_split_simulation": true,
    "seed": null
  }
}
```

### Fields

| Field | Type | Description | Default |
|-------|------|-------------|---------|
| `training_data_path` | string | NanoSim pre-trained model path | Required |
| `coverage` | integer | Target coverage depth | 50 |
| `min_read_length` | integer | Minimum read length (bp) | 1000 |
| `max_read_length` | integer | Maximum read length (bp) | 10000 |
| `correction_factor` | float | Coverage adjustment factor | 0.325 |
| `enable_split_simulation` | boolean | Diploid split-simulation mode | true |
| `seed` | integer | Random seed (null = random) | null |

### Diploid Split-Simulation

When `enable_split_simulation: true` and reference has 2 sequences:

1. Split diploid reference into haplotype1.fa and haplotype2.fa
2. Simulate each independently (coverage/2 each)
3. Merge reads from both haplotypes
4. Align merged reads to diploid reference

**Result:** Balanced allelic coverage (eliminates length-proportional bias).

---

## PacBio Parameters Section

### Purpose

Parameters for PacBio HiFi read simulation.

### Structure

```json
{
  "pacbio_params": {
    "model_type": "QSHMM",
    "model_file": "/path/to/pbsim3/models/QSHMM-RSII.model",
    "coverage": 30,
    "min_pass": 3,
    "max_pass": 15,
    "seed": null
  }
}
```

### Fields

| Field | Type | Description | Default |
|-------|------|-------------|---------|
| `model_type` | string | "QSHMM" or "ERRHMM" | "QSHMM" |
| `model_file` | string | pbsim3 model file path | Required |
| `coverage` | integer | Target coverage depth | 30 |
| `min_pass` | integer | Minimum CCS passes | 3 |
| `max_pass` | integer | Maximum CCS passes | 15 |
| `seed` | integer | Random seed (null = random) | null |

---

## Complete Example

```json
{
  "repeats": {
    "1": "AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG",
    "2": "AGTATGACCAGCAGCGTACTCTCCAGCCACAGCCCCGGTTCAGGCTCCTCCACCACTCAG",
    "X": "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA",
    "A": "GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA",
    "B": "GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA",
    "6": "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCGGGCCCCGGGCTCCACCCCGGCCCCG",
    "6p": "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCCCGGCCCCG",
    "7": "GGCTCCACCGCCCCCCCAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG",
    "8": "GGCTCCACCGCCCCCCCAGCCCATGGTGTCACCTCGGCCCCGGACAACAGGCCCGCCTTG",
    "9": "GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA"
  },

  "constants": {
    "hg38": {
      "left": "AGCAGGCAGTGCGGGGCCGCTGCTGCTG...",
      "right": "TGCTGCTGCTGCGGGGCCGCTGCTGCTG...",
      "vntr_start": 1500,
      "vntr_end": 8500
    }
  },

  "probabilities": {
    "1": {"2": 0.3, "X": 0.7},
    "2": {"1": 0.2, "A": 0.8},
    "X": {"A": 0.5, "B": 0.5},
    "A": {"B": 0.6, "X": 0.4},
    "B": {"A": 0.5, "X": 0.5}
  },

  "length_model": {
    "distribution_type": "normal",
    "mean_repeats": 63.3,
    "median_repeats": 70,
    "min_repeats": 42,
    "max_repeats": 85
  },

  "mutations": {
    "dupC": {
      "allowed_repeats": ["X", "A", "B"],
      "strict_mode": false,
      "changes": [
        {
          "operation": "insert",
          "sequence": "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA",
          "position_offset": 0
        }
      ]
    }
  },

  "tools": {
    "reseq": "/usr/bin/reseq",
    "bwa": "/usr/bin/bwa",
    "samtools": "/usr/bin/samtools"
  },

  "read_simulation": {
    "simulator": "illumina",
    "read_length": 150,
    "fragment_size": 350,
    "fragment_sd": 50,
    "coverage": 100,
    "threads": 4,
    "reference_genome": "/path/to/hg38.fa",
    "seed": 42
  },

  "nanosim_params": {
    "training_data_path": "/path/to/nanosim_training",
    "coverage": 50,
    "min_read_length": 1500,
    "max_read_length": 5000,
    "correction_factor": 0.325,
    "enable_split_simulation": true,
    "seed": 42
  },

  "pacbio_params": {
    "model_type": "QSHMM",
    "model_file": "/path/to/pbsim3/QSHMM-RSII.model",
    "coverage": 30,
    "min_pass": 3,
    "max_pass": 15,
    "seed": 42
  }
}
```

---

## Validation

MucOneUp validates configuration on load:

**Checked:**

- All repeat symbols referenced in probabilities exist in repeats
- Terminal block symbols (6/6p, 7, 8, 9) defined
- Probability values are valid (0.0 to 1.0)
- Length model parameters are positive integers
- Mutation sequences contain valid DNA (ACGT)
- Tool paths exist (if specified)

**Errors cause simulation to fail immediately.**

---

## Best Practices

!!! tip "Version Control Configuration"
    Commit `config.json` to version control with your simulation scripts for reproducibility.

!!! tip "Comment Your Mutations"
    JSON doesn't support comments, but you can use descriptive mutation names:
    ```json
    {
      "mutations": {
        "dupC_gastric_cancer_pmid12345": { ... }
      }
    }
    ```

!!! tip "Test Configuration"
    Validate configuration before large-scale simulations:
    ```bash
    muconeup --config config.json simulate --fixed-lengths 20 --out-base test
    ```

!!! warning "Platform-Specific Paths"
    Tool paths differ across systems. Use environment variables or separate configs for different machines:
    ```bash
    # Linux
    muconeup --config config_linux.json simulate ...

    # macOS
    muconeup --config config_macos.json simulate ...
    ```

---

## Next Steps

- **[Simulation Guide](../guides/simulation.md)** - Use your configuration for simulations
- **[Mutation Guide](../guides/mutations.md)** - Define custom mutations
- **[VNTR Statistics](../workflows/population-diversity.md)** - Analyze real data to inform probabilities

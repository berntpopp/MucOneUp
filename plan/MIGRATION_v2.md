# Migration Guide: v1.x → v2.0 (Click CLI)

## Overview

MucOneUp v2.0 introduces a new Click-based CLI with improved command structure following Unix philosophy. This guide provides side-by-side command comparisons to help you migrate from v1.x (argparse) to v2.0 (Click).

**Key Changes:**
- Single `muconeup` command → Four specialized commands: `simulate`, `reads`, `analyze`, `pipeline`
- Each command does ONE thing well (Unix philosophy)
- Cleaner command composition and chaining
- All v1.x functionality is preserved (100% feature parity)

## Command Architecture Changes

### v1.x (Argparse) - Single Command
```bash
muconeup [OPTIONS]
```
All functionality in one monolithic command with 25+ flags.

### v2.0 (Click) - Four Commands
```bash
muconeup simulate [OPTIONS]  # Generate VNTR haplotypes
muconeup reads [OPTIONS]     # Simulate sequencing reads
muconeup analyze [OPTIONS]   # Analyze sequences (ORFs, stats)
muconeup pipeline [OPTIONS]  # End-to-end automation
```

## Migration Examples

### Basic Haplotype Generation

**v1.x:**
```bash
muconeup --config config.json --out-base sample --out-dir output/
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample --out-dir output/
```

---

### Fixed Length Simulation

**v1.x:**
```bash
muconeup --config config.json --out-base sample --fixed-lengths 60
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample --fixed-lengths 60
```

---

### Length Range Series

**v1.x:**
```bash
muconeup --config config.json --out-base sample --fixed-lengths 20-40 --simulate-series
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample \
  --fixed-lengths 20-40 --simulate-series
```

---

### Mutation Application

**v1.x:**
```bash
muconeup --config config.json --out-base sample \
  --mutation-name dupC --mutation-targets 1,25 2,30
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample \
  --mutation-name dupC --mutation-targets 1,25 2,30
```

---

### Dual Simulation (Normal + Mutated)

**v1.x:**
```bash
muconeup --config config.json --out-base sample --mutation-name normal,dupC
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample --mutation-name normal,dupC
```

---

### Structure File Input

**v1.x:**
```bash
muconeup --config config.json --out-base sample --input-structure haplotypes.txt
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample --input-structure haplotypes.txt
```

---

### Random SNPs

**v1.x:**
```bash
muconeup --config config.json --out-base sample \
  --random-snps --random-snp-density 1.0
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample \
  --random-snps --random-snp-density 1.0
```

---

### SNPs from File

**v1.x:**
```bash
muconeup --config config.json --out-base sample --snp-file snps.tsv
```

**v2.0:**
```bash
muconeup simulate --config config.json --out-base sample --snp-file snps.tsv
```

---

### Read Simulation (Illumina)

**v1.x:**
```bash
muconeup --config config.json --out-base sample --simulate-reads
```

**v2.0 (Option 1 - Separate Commands):**
```bash
# Step 1: Generate haplotypes
muconeup simulate --config config.json --out-base sample --out-dir output/

# Step 2: Simulate reads
muconeup reads --config config.json --input output/sample.fa --out-dir output/
```

**v2.0 (Option 2 - Pipeline):**
```bash
muconeup pipeline illumina --config config.json --out-base sample --out-dir output/
```

---

### ONT Read Simulation

**v1.x:**
```bash
muconeup --config config.json --out-base sample --simulate-ont
```

**v2.0 (Option 1 - Separate Commands):**
```bash
# Step 1: Generate haplotypes
muconeup simulate --config config.json --out-base sample --out-dir output/

# Step 2: Simulate ONT reads
muconeup reads --config config.json --input output/sample.fa \
  --platform ont --out-dir output/
```

**v2.0 (Option 2 - Pipeline):**
```bash
muconeup pipeline ont --config config.json --out-base sample --out-dir output/
```

---

### ORF Prediction

**v1.x:**
```bash
muconeup --config config.json --out-base sample --output-orfs --orf-min-aa 100
```

**v2.0 (Option 1 - Separate Commands):**
```bash
# Step 1: Generate haplotypes
muconeup simulate --config config.json --out-base sample --out-dir output/

# Step 2: Predict ORFs
muconeup analyze orfs --input output/sample.fa --min-aa 100 --out-dir output/
```

**v2.0 (Option 2 - Pipeline):**
```bash
muconeup pipeline illumina --config config.json --out-base sample \
  --output-orfs --orf-min-aa 100 --out-dir output/
```

---

### Complex Workflow Example

**v1.x:**
```bash
muconeup --config config.json --out-base sample --out-dir output/ \
  --mutation-name dupC --mutation-targets 1,25 \
  --random-snps --random-snp-density 0.5 \
  --simulate-reads \
  --output-orfs --orf-min-aa 100 \
  --coverage 30 --threads 8 \
  --log-level INFO
```

**v2.0 (Option 1 - Step-by-Step):**
```bash
# Step 1: Generate mutated haplotypes with SNPs
muconeup simulate --config config.json --out-base sample --out-dir output/ \
  --mutation-name dupC --mutation-targets 1,25 \
  --random-snps --random-snp-density 0.5 \
  --log-level INFO

# Step 2: Simulate reads
muconeup reads --config config.json --input output/sample.fa --out-dir output/ \
  --coverage 30 --threads 8 --log-level INFO

# Step 3: Analyze ORFs
muconeup analyze orfs --input output/sample.fa --min-aa 100 --out-dir output/ \
  --log-level INFO
```

**v2.0 (Option 2 - Pipeline):**
```bash
muconeup pipeline illumina --config config.json --out-base sample --out-dir output/ \
  --mutation-name dupC --mutation-targets 1,25 \
  --random-snps --random-snp-density 0.5 \
  --output-orfs --orf-min-aa 100 \
  --coverage 30 --threads 8 \
  --log-level INFO
```

---

## Complete Option Mapping

### Configuration Options
| v1.x Flag | v2.0 Command | v2.0 Flag | Notes |
|-----------|--------------|-----------|-------|
| `--config` | `simulate` / `reads` / `pipeline` | `--config` | Unchanged |
| `--out-base` | `simulate` / `pipeline` | `--out-base` | Unchanged |
| `--out-dir` | All commands | `--out-dir` | Unchanged |
| `--log-level` | All commands | `--log-level` | Unchanged |

### Haplotype Generation Options
| v1.x Flag | v2.0 Command | v2.0 Flag | Notes |
|-----------|--------------|-----------|-------|
| `--fixed-lengths` | `simulate` | `--fixed-lengths` | Unchanged |
| `--simulate-series` | `simulate` | `--simulate-series` | Unchanged |
| `--num-iterations` | `simulate` | `--num-iterations` | Unchanged |
| `--input-structure` | `simulate` | `--input-structure` | Unchanged |
| `--input-chains` | `simulate` | `--input-chains` | Unchanged (alias) |

### Mutation Options
| v1.x Flag | v2.0 Command | v2.0 Flag | Notes |
|-----------|--------------|-----------|-------|
| `--mutation-name` | `simulate` / `pipeline` | `--mutation-name` | Unchanged |
| `--mutation-targets` | `simulate` / `pipeline` | `--mutation-targets` | Unchanged |

### SNP Options
| v1.x Flag | v2.0 Command | v2.0 Flag | Notes |
|-----------|--------------|-----------|-------|
| `--snp-file` | `simulate` / `pipeline` | `--snp-file` | Unchanged |
| `--random-snps` | `simulate` / `pipeline` | `--random-snps` | Unchanged |
| `--random-snp-density` | `simulate` / `pipeline` | `--random-snp-density` | Unchanged |

### Read Simulation Options
| v1.x Flag | v2.0 Command | v2.0 Flag | Notes |
|-----------|--------------|-----------|-------|
| `--simulate-reads` | `reads` / `pipeline` | N/A | Command implies action |
| `--simulate-ont` | `reads --platform ont` | `--platform ont` | Now explicit flag |
| `--coverage` | `reads` / `pipeline` | `--coverage` | Unchanged |
| `--threads` | `reads` / `pipeline` | `--threads` | Unchanged |
| `--downsample` | `reads` / `pipeline` | `--downsample` | Unchanged |

### Analysis Options
| v1.x Flag | v2.0 Command | v2.0 Flag | Notes |
|-----------|--------------|-----------|-------|
| `--output-orfs` | `analyze orfs` / `pipeline` | N/A | Command implies action |
| `--orf-min-aa` | `analyze orfs` / `pipeline` | `--min-aa` | Renamed for clarity |
| `--orf-all-frames` | `analyze orfs` / `pipeline` | `--all-frames` | Prefix removed |

## When to Use Each Command

### Use `simulate` When:
- ✅ You only need FASTA haplotypes (no reads)
- ✅ You want to inspect sequences before read simulation
- ✅ You're generating training data
- ✅ You need structure files for downstream analysis

### Use `reads` When:
- ✅ You already have FASTA input (from `simulate` or elsewhere)
- ✅ You want to simulate reads from existing sequences
- ✅ You need fine control over read simulation parameters

### Use `analyze` When:
- ✅ You want to analyze existing FASTA files
- ✅ You need ORF prediction only
- ✅ You want statistics without full simulation

### Use `pipeline` When:
- ✅ You want end-to-end automation (haplotypes → reads → analysis)
- ✅ You're running production workflows
- ✅ You need consistent parameter application across all steps
- ✅ You want single-command convenience

## Command Chaining Patterns

### Pattern 1: Generate Once, Simulate Multiple Times
```bash
# Generate haplotypes once
muconeup simulate --config config.json --out-base sample --out-dir output/

# Simulate different coverage depths
muconeup reads --config config.json --input output/sample.fa \
  --coverage 10 --out-dir output/low/
muconeup reads --config config.json --input output/sample.fa \
  --coverage 50 --out-dir output/high/
```

### Pattern 2: Compare Platforms
```bash
# Generate haplotypes once
muconeup simulate --config config.json --out-base sample --out-dir output/

# Simulate both Illumina and ONT
muconeup reads --config config.json --input output/sample.fa \
  --platform illumina --out-dir output/illumina/
muconeup reads --config config.json --input output/sample.fa \
  --platform ont --out-dir output/ont/
```

### Pattern 3: Analyze Existing Data
```bash
# Analyze pre-existing FASTA
muconeup analyze orfs --input existing_sample.fa --out-dir analysis/
```

### Pattern 4: Full Automation
```bash
# Single command for complete workflow
muconeup pipeline illumina --config config.json --out-base sample \
  --out-dir output/ --output-orfs --coverage 30
```

## Breaking Changes

### ⚠️ None! (100% Backward Compatible)

All v1.x functionality is preserved in v2.0. The changes are architectural:
- Old command structure → New command structure
- Single command → Multiple specialized commands
- Flags are identical, just organized under appropriate commands

### Command Line Length
v2.0 commands may be longer for complex workflows when using separate commands. Use `pipeline` for convenience or shell scripts for automation.

## Shell Script Migration

### v1.x Script
```bash
#!/bin/bash
for length in 20 40 60 80; do
  muconeup --config config.json --out-base sample_${length} \
    --fixed-lengths $length --simulate-reads --coverage 30
done
```

### v2.0 Script (Option 1: Pipeline)
```bash
#!/bin/bash
for length in 20 40 60 80; do
  muconeup pipeline illumina --config config.json --out-base sample_${length} \
    --fixed-lengths $length --coverage 30 --out-dir output/
done
```

### v2.0 Script (Option 2: Separate Steps)
```bash
#!/bin/bash
for length in 20 40 60 80; do
  # Generate haplotypes
  muconeup simulate --config config.json --out-base sample_${length} \
    --fixed-lengths $length --out-dir output/

  # Simulate reads
  muconeup reads --config config.json --input output/sample_${length}.fa \
    --coverage 30 --out-dir output/
done
```

## Help and Documentation

### Get Help
```bash
# Top-level help
muconeup --help

# Command-specific help
muconeup simulate --help
muconeup reads --help
muconeup analyze --help
muconeup pipeline --help
```

### Version Check
```bash
# Both v1.x and v2.0
muconeup --version
```

## Troubleshooting

### "Command not found: muconeup simulate"
**Cause:** You're running v1.x, not v2.0.

**Solution:** Upgrade to v2.0:
```bash
pip install --upgrade muconeup
```

### "Too many arguments"
**Cause:** You used v1.x syntax with v2.0 installation.

**Solution:** Add command name after `muconeup`:
```bash
# Old (v1.x)
muconeup --config config.json --out-base sample

# New (v2.0)
muconeup simulate --config config.json --out-base sample
```

### "No such option: --simulate-reads"
**Cause:** You used v1.x flag with v2.0 `simulate` command.

**Solution:** Either:
1. Use `reads` command separately
2. Use `pipeline` command for end-to-end

```bash
# Option 1: Separate commands
muconeup simulate --config config.json --out-base sample --out-dir output/
muconeup reads --config config.json --input output/sample.fa --out-dir output/

# Option 2: Pipeline
muconeup pipeline illumina --config config.json --out-base sample --out-dir output/
```

## Benefits of v2.0

### 1. Unix Philosophy
Each command does ONE thing well:
- `simulate` → Generate sequences
- `reads` → Simulate reads
- `analyze` → Analyze data
- `pipeline` → Automate workflows

### 2. Composability
Commands can be chained flexibly:
```bash
muconeup simulate ... && muconeup reads ... && muconeup analyze ...
```

### 3. Clarity
Clear command names make intent obvious:
```bash
muconeup reads --platform ont  # Clear: simulating ONT reads
```
vs.
```bash
muconeup --simulate-ont  # Less clear: what is being done?
```

### 4. Extensibility
New analysis commands can be added cleanly:
```bash
muconeup analyze variants ...  # Future command
muconeup analyze coverage ...  # Future command
```

### 5. Help Clarity
Each command has focused help:
```bash
muconeup simulate --help  # Only haplotype generation options
muconeup reads --help     # Only read simulation options
```
vs.
```bash
muconeup --help  # 25+ options mixed together (v1.x)
```

## Support

- **Documentation:** See README.md for detailed examples
- **Issues:** https://github.com/YOUR_ORG/MucOneUp/issues
- **Questions:** Open a GitHub discussion

## Version Compatibility

- **v1.x:** Monolithic argparse CLI (deprecated)
- **v2.0:** Click-based CLI with command structure (current)
- **Forward:** All future features will use v2.0 architecture

**Recommendation:** Migrate to v2.0 for long-term support and new features.

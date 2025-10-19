# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**MucOneUp** is a Python tool for simulating MUC1 VNTR (Variable Number Tandem Repeat) diploid references with customizable mutations and read simulation. The tool generates realistic genomic sequences with configurable repeat structures, applies targeted mutations, and can simulate both Illumina and Oxford Nanopore sequencing reads.

## Common Commands

### Installation
```bash
pip install .
```

### Running Tests
```bash
# Run all tests
python -m pytest tests/

# Run specific test file
python -m pytest tests/test_simulate.py

# Run with verbose output (pytest)
python -m pytest -v tests/

# Run with MucOneUp verbose flag
muconeup --verbose --config config.json simulate --help
muconeup -v --config config.json simulate --help  # Short form
```

### Basic Simulation
```bash
# Generate diploid haplotypes with random VNTR lengths
muconeup --config config.json simulate --out-base output_name --out-dir output/

# Generate with fixed VNTR lengths
muconeup --config config.json simulate --out-base output_name --fixed-lengths 60

# Generate series across length range (with progress bar)
muconeup --config config.json simulate --out-base output_name --fixed-lengths 20-40 --simulate-series 5

# With verbose output
muconeup --verbose --config config.json simulate --out-base output_name --fixed-lengths 20-40

# Apply mutations to specific positions
muconeup --config config.json simulate --out-base output_name --mutation-name dupC --mutation-targets 1,25 2,30

# Dual simulation (normal + mutated)
muconeup --config config.json simulate --out-base output_name --mutation-name normal,dupC
```

### Progress Indicators
When using `--simulate-series` with multiple iterations, a progress bar will appear:
```
Simulating 21 iterations  [################-----------]   61%  00:02:15
```
This helps track long-running simulations.

### Read Simulation
```bash
# Simulate Illumina reads
muconeup --config config.json --out-base output_name --simulate-reads

# With ORF prediction
muconeup --config config.json --out-base output_name --output-orfs --orf-min-aa 100
```

### Deterministic Read Simulation

Generate reproducible reads using the `--seed` parameter. This ensures identical read output across multiple runs, enabling reproducible research, fair benchmarking, and consistent debugging.

```bash
# Illumina reads with seed (reproducible)
muconeup --config config.json reads illumina sample.fa --seed 42

# ONT reads with seed (reproducible)
muconeup --config config.json reads ont sample.fa --seed 42

# Full pipeline with seed (haplotypes + reads)
muconeup --config config.json simulate --seed 42 --out-base sample
muconeup --config config.json reads illumina sample.001.simulated.fa --seed 42

# Verify reproducibility
muconeup --config config.json reads illumina test.fa --seed 42 --out-base run1
muconeup --config config.json reads illumina test.fa --seed 42 --out-base run2
diff run1_R1.fastq.gz run2_R1.fastq.gz  # Should be identical
```

**Important**: Same seed guarantees identical output ONLY when:
- Using identical input files
- Running on same platform/architecture
- Using same tool versions (NanoSim, Python, reseq)

**Configuration Example**:
```json
{
  "read_simulation": {
    "simulator": "illumina",
    "seed": 42
  },
  "nanosim_params": {
    "training_data_path": "/path/to/model",
    "coverage": 30,
    "seed": 42
  }
}
```

### SNP Integration
```bash
# Generate random SNPs
muconeup --config config.json --out-base output_name --random-snps --random-snp-density 1.0

# Apply SNPs from file
muconeup --config config.json --out-base output_name --snp-file snps.tsv
```

### Setting Up Conda Environments
```bash
# For Illumina read simulation (w-Wessim2)
mamba env create -f conda/env_wessim.yaml

# For ONT read simulation (NanoSim)
mamba env create -f conda/env_nanosim.yml
```

## Architecture

### Core Simulation Pipeline

1. **Configuration Loading** (`config.py`): Validates JSON configuration against schema containing repeat definitions, probability transitions, mutation definitions, and tool paths.

2. **VNTR Generation** (`simulate.py`):
   - Builds haplotype chains by sampling repeats according to probability distributions
   - Forces canonical terminal block (6/6p → 7 → 8 → 9)
   - Supports both random length sampling and fixed-length generation
   - Can generate from predefined structure files

3. **Mutation Application** (`mutate.py`):
   - Applies insertions, deletions, replacements, or delete_insert operations
   - Supports strict mode to prevent auto-conversion of non-allowed repeats
   - Tracks mutated positions with "m" suffix in structure files
   - Records mutated VNTR unit sequences separately

4. **SNP Integration** (`snp_integrator.py`):
   - Parses TSV files with haplotype-specific SNPs (1-based haplotype, 0-based position)
   - Generates random SNPs with configurable density
   - Validates reference bases before applying (skippable in dual mutation mode)
   - Tracks successfully applied SNPs for reporting

5. **Sequence Assembly**: Concatenates left constant → repeat chain → right constant, supporting both hg19 and hg38 assemblies

### Read Simulation Pipelines

#### Illumina Pipeline (`read_simulator/pipeline.py`)
Uses a port of w-Wessim2 with these steps:
1. Replace Ns using reseq
2. Generate systematic errors with reseq illuminaPE
3. Convert to 2bit format with faToTwoBit
4. Extract subset reference from sample BAM
5. Align with pblat
6. Simulate fragments (ported w-Wessim2 logic in `fragment_simulation.py`)
7. Create paired reads with reseq seqToIllumina
8. Split interleaved FASTQ
9. Align to human reference with BWA MEM
10. Optional coverage-based downsampling

#### ONT Pipeline (`read_simulator/ont_pipeline.py`)
Uses NanoSim for Oxford Nanopore long reads with automatic diploid split-simulation support.

**Standard Mode (Haploid/Single Sequence)**:
1. Run NanoSim simulation with pre-trained models
2. Apply coverage correction factor (default: 0.325)
3. Align reads with minimap2
4. Create indexed BAM output

**Diploid Split-Simulation Mode (2 Sequences)** - NEW in v0.14.0:
Automatically activated when reference contains exactly 2 sequences. Eliminates allelic bias in diploid references.

Workflow:
1. **Diploid Detection**: Automatically identifies 2-sequence FASTA files
2. **Haplotype Extraction**: Splits diploid reference into separate haplotype files
3. **Coverage Calculation**: Divides target coverage by 2, applies correction factor
4. **Independent Simulation**: Simulates each haplotype separately with unique seeds (seed, seed+1)
5. **FASTQ Merging**: Combines reads from both haplotypes
6. **Alignment**: Aligns merged reads to reference
7. **BAM Output**: Creates indexed BAM file

Features:
- Automatic mode detection (no user intervention)
- Eliminates length-proportional sampling bias
- Maintains reproducibility with seed support
- Backward compatible with haploid references

**Read Length Control** (recommended for diploid):
```json
"nanosim_params": {
  "min_read_length": 1500,
  "max_read_length": 5000,
  "correction_factor": 0.325,
  "enable_split_simulation": true
}
```

### Key Modules

- **`cli.py`**: Argument parsing, simulation orchestration, dual simulation mode, series generation
- **`distribution.py`**: Samples target VNTR length from normal/uniform distributions
- **`probabilities.py`**: Weighted random selection for repeat transitions
- **`fasta_writer.py`**: FASTA output with per-haplotype mutation annotations
- **`io.py`**: Structure file parsing with embedded mutation information
- **`translate.py`**: ORF prediction using orfipy
- **`toxic_protein_detector.py`**: Scans ORFs for toxic features based on repeat structure and amino acid composition
- **`simulation_statistics.py`**: Generates comprehensive JSON reports with runtime, haplotype metrics, and mutation details
- **`analysis/vntr_statistics.py`**: Analyzes VNTR structures from CSV/TSV files, computes statistics (min/max/mean/median repeats), and builds transition probability matrices

### Configuration Structure

The `config.json` file contains:
- **repeats**: Dictionary mapping repeat symbols (1, 2, X, A, B, etc.) to DNA sequences
- **constants**: Left/right flanking sequences for hg19 and hg38, plus VNTR region coordinates
- **probabilities**: State transition probabilities for repeat chain generation
- **length_model**: Distribution parameters (normal/uniform) with min/max/mean/median repeats
- **mutations**: Named mutation definitions with:
  - `allowed_repeats`: Valid repeat symbols for this mutation
  - `strict_mode`: Boolean to enforce allowed_repeats (prevents auto-conversion)
  - `changes`: List of operations (insert/delete/replace/delete_insert) with positions and sequences
- **tools**: Command paths for external tools (reseq, bwa, samtools, etc.)
- **read_simulation**: Parameters for Illumina pipeline (fragment size, coverage, threads, etc.)
- **nanosim_params**: Parameters for ONT pipeline (training model path, coverage, read lengths)

### Mutation Strict Mode

By default, if a mutation target has a repeat not in `allowed_repeats`, it's auto-converted with a warning. Setting `"strict_mode": true` in a mutation definition causes an error instead. Random targets always respect `allowed_repeats` in both modes.

### Structure Files

Structure files can contain mutation information in header comments:
```
# Mutation Applied: dupC (Targets: [(1, 25)])
haplotype_1 1-2-3-4-5-C-X-B-Xm-X-A-6-7-8-9
haplotype_2 1-2-3-4-5-C-X-A-B-X-6p-7-8-9
```

The "m" suffix marks mutated positions. Use `--input-structure` to generate from predefined chains.

### Dual Simulation Mode

When `--mutation-name normal,mutationName` is provided, two complete simulation runs occur:
- `*.normal.fa` and `*.normal.simulation_stats.json`
- `*.mut.fa` and `*.mut.simulation_stats.json`

In dual mode, SNP integration uses `skip_reference_check=True` for mutated sequences.

### Reference Assembly Support

The tool supports both hg19 and hg38 assemblies. Set `"reference_assembly": "hg38"` or `"hg19"` in config. Constants and VNTR regions are assembly-specific.

## Development Notes

- **Entry point**: `muconeup` command maps to `muc_one_up.cli:main` (defined in setup.cfg)
- **Version management**: Single source in `muc_one_up/version.py`, imported by `__init__.py` and `cli.py`
- **File naming**: Simulation outputs use numbered iterations (`.001`, `.002`, etc.) for series mode
- **External tools**: All pipeline wrappers in `read_simulator/wrappers/` have timeout handling
- **Intermediate files**: Pipeline uses underscore prefix (`_`) for temporary files
- **Logging**: Configurable via `--log-level` (DEBUG/INFO/WARNING/ERROR/CRITICAL/NONE)
- **Example data**: Located in `data/examples/`, includes `vntr_database.tsv` with real-world VNTR structures from published research

## Testing Strategy

Tests cover:
- Configuration validation and loading (`test_config.py`)
- Probability-based repeat selection (`test_probabilities.py`)
- Haplotype simulation and chain building (`test_simulate.py`)
- Mutation application and validation (`test_mutate.py`)
- ORF prediction and translation (`test_translate.py`)
- VNTR structure analysis (`test_vntr_statistics.py`)
- CLI integration for vntr-stats command (`test_click_vntr_stats.py`)

When adding features, ensure corresponding tests validate both success and error cases.

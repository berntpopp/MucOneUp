# MucOneUp

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/MucOneUp/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneUp/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **MUC1 VNTR simulation and analysis toolkit for genomics research**

---

## Overview

MucOneUp is a **Python toolkit for simulating realistic MUC1 Variable Number Tandem Repeat (VNTR) sequences** with customizable mutations and sequencing read generation. Designed for genomics researchers studying MUC1 gene variation, it enables reproducible generation of diploid haplotypes with targeted mutations, SNP integration, and comprehensive downstream analysis including ORF prediction and toxic protein detection.

**Perfect for:** Benchmarking variant callers, testing mutation detection pipelines, generating synthetic training data, and exploring MUC1 VNTR structural diversity.

---

## Key Features

- 🧬 **Realistic VNTR Simulation** - Probability-based repeat transitions with canonical terminal blocks (6/6p → 7 → 8 → 9)
- 🔬 **Flexible Mutation Engine** - Insert, delete, replace, or delete-insert operations with strict mode validation
- 📊 **Multi-Platform Read Simulation** - Illumina (w-Wessim2) and Oxford Nanopore (NanoSim) integration
- 🧪 **ORF Prediction & Toxic Detection** - Automated open reading frame analysis with toxicity scoring
- 🧮 **SNP Integration** - Random or predefined SNP application with haplotype-specific variants
- 🔄 **Batch Processing** - Unix-style composable commands (`simulate` → `analyze` → `reads`)
- 📈 **Comprehensive Statistics** - JSON reports with per-haplotype metrics and mutation tracking

---

## Core Design Principles

MucOneUp follows modern software engineering best practices:

- **Unix Philosophy** - Each command does ONE thing well; compose for complex workflows
- **SOLID Architecture** - Modular, testable components with clear separation of concerns
- **Type Safety** - Comprehensive type hints with mypy validation
- **DRY (Don't Repeat Yourself)** - Centralized utilities, zero code duplication
- **KISS (Keep It Simple)** - Minimal configuration, sensible defaults
- **Quality First** - 100% Google-style docstrings, 77% test coverage, zero linting errors

**Built with:** Python 3.10+, Click CLI, pytest, ruff, mypy, pre-commit hooks

---

## Installation

### Quick Install (Users)

```bash
pip install .
```

### Development Setup

Modern Python tooling with **uv**, **ruff**, **mypy**, and automated **pre-commit hooks**:

```bash
# Install uv (fast package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Setup environment and dependencies
make init

# Verify installation
make check
```

**Development Commands:**

| Command | Action |
|---------|--------|
| `make init` | Setup complete dev environment |
| `make test` | Run tests with coverage (568 tests) |
| `make lint` | Check code quality (ruff + mypy) |
| `make format` | Auto-format code |
| `make check` | Run all quality checks |

### Additional Setup

- 📚 **Reference Files** - Human genome (hg19/hg38), reseq models - See `helpers/install_references.py`
- 🐚 **Shell Completion** - Bash/Zsh/Fish autocomplete - Run `_MUCONEUP_COMPLETE=bash_source muconeup`
- 🧪 **Read Simulation** - Conda/Mamba environments - See `conda/env_wessim.yaml` and `conda/env_nanosim.yml`

**System Requirements:** Python 3.10+, 4GB RAM minimum, 50GB disk for reference files

---

## Quick Start

Generate a diploid haplotype with mutation and analyze the results:

```bash
# 1. Generate diploid haplotypes with mutation
muconeup --config config.json simulate \
  --out-base muc1_example \
  --out-dir output/ \
  --mutation-name dupC \
  --mutation-targets 1,25 \
  --output-structure

# 2. Predict ORFs and detect toxic proteins
muconeup --config config.json analyze orfs \
  output/muc1_example.001.simulated.fa \
  --out-base muc1_orfs \
  --orf-min-aa 100

# 3. Simulate Illumina reads
muconeup --config config.json reads illumina \
  output/muc1_example.001.simulated.fa \
  --out-base muc1_reads \
  --coverage 100
```

**Output Files:**
```
output/
├── muc1_example.001.simulated.fa          # Diploid haplotype sequences
├── muc1_example.001.vntr_structure.txt    # Repeat chain structure
├── muc1_example.001.simulation_stats.json # Comprehensive metrics
├── muc1_orfs.pep.fa                       # Predicted ORF peptides
├── muc1_orfs.orf_stats.txt                # Toxic protein detection
└── muc1_reads.illumina.bam                # Simulated reads (aligned)
```

**What This Example Shows:**
- ✅ Haplotype generation with fixed mutation
- ✅ Command composition (Unix philosophy)
- ✅ Complete workflow from simulation → analysis → reads
- ✅ Output file organization

**Next Steps:** See full documentation for advanced workflows, batch processing, and configuration options.

---

## Documentation

### User Guides

- **[CLI Reference](plan/MIGRATION_v2.md)** - All commands and options (migration guide until Sphinx docs complete)
- **[Configuration Guide](config.json)** - Example config with mutations, probabilities, and tool paths
- **[Installation Details](helpers/install_references.py)** - Reference file setup helper script
- **[Development Plan](plan/README.md)** - Complete modernization history and roadmap

### CLI Commands

**Core Commands:**
- `muconeup simulate` - Generate haplotypes ONLY (core functionality)
- `muconeup reads {illumina|ont}` - Simulate sequencing reads from ANY FASTA
- `muconeup analyze {orfs|stats|vntr-stats}` - Analyze ANY FASTA file

**Batch Processing:**
```bash
# Multiple files at once
muconeup reads illumina *.simulated.fa --coverage 30

# Or use xargs/GNU parallel
ls *.simulated.fa | parallel muconeup reads illumina {}
```

### Advanced Topics

- **Dual Simulation Mode** - Generate normal + mutated pairs with `--mutation-name normal,dupC`
- **Structure Files** - Reproducible mutations from predefined repeat chains
- **SNP Integration** - Random (`--random-snps`) or file-based (`--snp-input-file`) variants
- **Series Generation** - Sweep length ranges with `--fixed-lengths 20-40 --simulate-series 5`
- **Mutation Strict Mode** - Enforce allowed repeats with `"strict_mode": true` in config
- **VNTR Statistics** - Analyze repeat databases with `muconeup analyze vntr-stats`

See `plan/completed/` for detailed implementation documentation.

---

## Development

### Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Code style guidelines (ruff, mypy, Google-style docstrings)
- Testing requirements (pytest, 77% coverage minimum)
- Pull request process

### Development Workflow

```bash
make init        # Setup dev environment
make test        # Run test suite (568 tests)
make lint        # Check code quality (ruff + mypy)
make format      # Auto-format code
make check       # Run all quality checks
```

### Quality Standards

- ✅ **100% Google-style docstring coverage** (v0.12.0 milestone)
- ✅ **77% test coverage** - 568 tests passing across all modules
- ✅ **Zero linting errors** - ruff compliance enforced by pre-commit hooks
- ✅ **Zero type errors** - mypy strict mode passing
- ✅ **Automated quality gates** - pre-commit hooks + GitHub Actions CI/CD

### Project Structure

```
muc_one_up/
├── cli/                    # Click command-line interface
│   ├── commands/           # Command implementations (simulate, reads, analyze)
│   ├── orchestration.py    # Workflow orchestration
│   └── ...                 # CLI utilities (mutations, outputs, SNPs)
├── bioinformatics/         # DNA/FASTA/SNP validation
├── read_simulator/         # Illumina & ONT pipelines
│   ├── pipeline.py         # Illumina w-Wessim2 integration
│   ├── ont_pipeline.py     # Oxford Nanopore NanoSim integration
│   └── wrappers/           # External tool wrappers (DRY/SOLID)
├── analysis/               # VNTR statistics, ORF prediction
├── simulate.py             # Core VNTR haplotype simulation engine
├── mutate.py               # Mutation application (insert/delete/replace)
├── probabilities.py        # Weighted random repeat selection
├── distribution.py         # Target length sampling
└── type_defs.py            # Type aliases & Protocol definitions
```

**Architecture Highlights:**
- **Centralized command utilities** - `command_utils.py` eliminates 96% code duplication
- **DRY/SOLID wrappers** - All external tools use centralized, security-compliant builders
- **100% test coverage** on critical modules (config, probabilities, validation)

See [plan/README.md](plan/README.md) for complete refactoring history (9 completed phases).

---

## Citation

If you use MucOneUp in your research, please cite:

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  version = {0.13.0},
  url = {https://github.com/berntpopp/MucOneUp}
}
```

**Development Status:** Pre-release software under active development.

A manuscript describing MucOneUp is in preparation. For now, please cite using the software reference above along with the GitHub repository URL and version number used in your research.

**Community:** If you publish results using MucOneUp, please let us know by [opening an issue](https://github.com/berntpopp/MucOneUp/issues)! We'd love to showcase your work.

---

## License

This project is licensed under the **MIT License** - see [LICENSE](LICENSE) file for details.

**In short:** ✅ Commercial use ✅ Modification ✅ Distribution ✅ Private use

---

**Version:** 0.13.0 | **Status:** Active Development | **Maintained By:** [Bernt Popp](https://github.com/berntpopp)

# MucOneUp

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/MucOneUp/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneUp/actions)
[![Docker](https://img.shields.io/badge/docker-ghcr.io-blue.svg)](https://github.com/berntpopp/MucOneUp/pkgs/container/muconeup%2Fmuconeup)
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
- 📊 **Multi-Platform Read Simulation** - Illumina (w-Wessim2), Oxford Nanopore (NanoSim), and PacBio HiFi (pbsim3/CCS) integration
- 🧪 **ORF Prediction & Toxic Detection** - Automated open reading frame analysis with toxicity scoring
- 🧮 **SNP Integration** - Random or predefined SNP application with haplotype-specific variants
- 🔄 **Batch Processing** - Unix-style composable commands (`simulate` → `analyze` → `reads`)
- 📈 **Comprehensive Statistics** - JSON reports with per-haplotype metrics and mutation tracking

---

## Installation

### Quick Install (Users)

```bash
make install
```

### Development Setup

Modern Python tooling with **uv**, **ruff**, **mypy**, and automated **pre-commit hooks**:

```bash
make init    # Installs uv, dev dependencies, pre-commit hooks
make check   # Verify installation
```

**Common Commands:**

| Command | Action |
|---------|--------|
| `make install` | Install package for users |
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

## 🐳 Docker Installation

**Simplest option** - No Python environment needed, all tools included:

```bash
# Pull image from GitHub Container Registry
docker pull ghcr.io/berntpopp/muconeup/muconeup:latest

# Run
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  ghcr.io/berntpopp/muconeup/muconeup:latest \
  --config /app/config.json \
  simulate --out-base /data/sample
```

**What's included:** MucOneUp core + Illumina (w-Wessim2) + ONT (NanoSim) + PacBio (pbsim3/CCS) - all in one image.

### Docker Compose

```bash
# Create docker-compose.yml (or use the one in repo)
docker-compose run --rm muconeup --config /app/config.json simulate --help
```

### Build from Source

```bash
./docker/build.sh
```

**Advantages:**
- ✅ No conda/Python setup required
- ✅ All simulators pre-installed
- ✅ Reproducible environment
- ✅ Easy deployment on HPC/cloud

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
- `muconeup reads {illumina|ont|pacbio}` - Simulate sequencing reads from ANY FASTA
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

## Citation

If you use MucOneUp in your research, please cite:

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  note = {Software version number available at: https://github.com/berntpopp/MucOneUp/releases},
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

**Status:** Active Development | **Maintained By:** [Bernt Popp](https://github.com/berntpopp)

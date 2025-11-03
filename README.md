# MucOneUp

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/MucOneUp/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneUp/actions)
[![Documentation](https://img.shields.io/badge/docs-MkDocs%20Material-blue)](https://berntpopp.github.io/MucOneUp/)
[![Docker](https://img.shields.io/badge/docker-ghcr.io-blue.svg)](https://github.com/berntpopp/MucOneUp/pkgs/container/muconeup%2Fmuconeup)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MUC1 VNTR simulation and analysis toolkit for genomics research.

---

## Overview

MucOneUp generates realistic **MUC1 Variable Number Tandem Repeat (VNTR) sequences** with customizable mutations and platform-specific sequencing reads. Designed for benchmarking variant callers, testing mutation detection pipelines, and generating synthetic training data.

**Key Capabilities:**

- Diploid haplotype generation with probability-based repeat transitions
- Frameshift mutation simulation (dupC, delC, custom insertions/deletions)
- Multi-platform read simulation (Illumina, Oxford Nanopore, PacBio HiFi)
- ORF prediction with toxic protein detection for ADTKD-MUC1 analysis
- SNP integration and VNTR statistics analysis
- Reproducible workflows with seed-based generation

---

## Quick Start

```bash
# View help and version (no --config required)
muconeup -h        # Show help
muconeup -V        # Show version

# Generate diploid haplotypes with mutation
muconeup --config config.json simulate \
  --out-base sample \
  --mutation-name dupC \
  --mutation-targets 1,25

# Predict ORFs and detect toxic proteins
muconeup --config config.json analyze orfs \
  sample.001.simulated.fa \
  --out-base orfs

# Simulate Illumina reads
muconeup --config config.json reads illumina \
  sample.001.simulated.fa \
  --coverage 100
```

**Documentation:** https://berntpopp.github.io/MucOneUp/

---

## Installation

### Option 1: pip (Recommended for Users)

```bash
pip install git+https://github.com/berntpopp/MucOneUp.git
```

**Requirements:**
- Python 3.10+
- External tools for read simulation (optional): BWA, samtools, reseq

### Option 2: Docker (All-in-One)

```bash
# Pull image
docker pull ghcr.io/berntpopp/muconeup/muconeup:latest

# Run
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  ghcr.io/berntpopp/muconeup/muconeup:latest \
  --config /app/config.json \
  simulate --out-base /data/sample
```

**Includes:** MucOneUp + Illumina (w-Wessim2) + ONT (NanoSim) + PacBio (pbsim3) - pre-configured environment.

### Option 3: Development Setup

```bash
# Clone repository
git clone https://github.com/berntpopp/MucOneUp.git
cd MucOneUp

# Install with uv (modern Python package manager)
make init    # Install dev dependencies + pre-commit hooks
make test    # Run 568 tests
make check   # Verify installation
```

---

## Citation

If you use MucOneUp in your research:

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  url = {https://github.com/berntpopp/MucOneUp},
  note = {Version available at https://github.com/berntpopp/MucOneUp/releases}
}
```

**Status:** Pre-release software under active development. A manuscript is in preparation.

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

**Maintained by:** [Bernt Popp](https://github.com/berntpopp)

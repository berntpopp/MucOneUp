# Installation

MucOneUp supports multiple installation methods for different use cases.

---

## Prerequisites

**System Requirements:**

- **Python:** 3.10, 3.11, 3.12, or 3.13
- **RAM:** 4GB minimum (8GB recommended for read simulation)
- **Disk Space:** 50GB for reference files (hg19/hg38, reseq models)
- **Operating System:** Linux, macOS, Windows (via WSL2)

**Required for Read Simulation:**

- **Illumina:** reseq, bwa, samtools, UCSC tools (faToTwoBit, pblat)
- **Oxford Nanopore:** NanoSim, minimap2, samtools
- **PacBio HiFi:** pbsim3, ccs, minimap2, samtools

---

## Quick Install (Recommended)

### Standard Installation

Install MucOneUp and core dependencies:

```bash
# Clone repository
git clone https://github.com/berntpopp/MucOneUp.git
cd MucOneUp

# Install package
make install

# Verify installation
muconeup --version
```

**What this installs:**

- MucOneUp CLI tool
- Core Python dependencies (Click, orfipy, jsonschema)
- Entry point for `muconeup` command

**Does NOT include:** External tools for read simulation (install separately via conda).

---

## Docker Installation (Easiest)

Pre-built Docker images include **all dependencies** (MucOneUp + Illumina/ONT/PacBio tools).

### Pull from GitHub Container Registry

```bash
# Pull latest image
docker pull ghcr.io/berntpopp/muconeup/muconeup:latest

# Verify installation
docker run --rm ghcr.io/berntpopp/muconeup/muconeup:latest --version
```

### Run Simulation

```bash
# Create local directories
mkdir -p data output

# Copy config to current directory
cp config.json .

# Run simulation (mount volumes)
docker run --rm \
  -v $(pwd)/output:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  ghcr.io/berntpopp/muconeup/muconeup:latest \
  --config /app/config.json \
  simulate --out-base /data/sample
```

### Docker Compose

Create `docker-compose.yml`:

```yaml
version: '3.8'

services:
  muconeup:
    image: ghcr.io/berntpopp/muconeup/muconeup:latest
    volumes:
      - ./data:/data
      - ./config.json:/app/config.json:ro
    command: --config /app/config.json simulate --out-base /data/sample
```

Run:

```bash
docker-compose run --rm muconeup --config /app/config.json simulate --help
```

**Advantages:**

- ✅ No conda/Python environment setup
- ✅ All simulators pre-installed
- ✅ Reproducible environment
- ✅ Easy HPC/cloud deployment

---

## Development Setup

For contributors and developers:

### Modern Python Tooling

MucOneUp uses **uv** (fast Python package installer), **ruff** (linter/formatter), and **mypy** (type checker).

```bash
# Install development environment
make init

# What this does:
# 1. Installs uv package manager
# 2. Installs development dependencies
# 3. Sets up pre-commit hooks (ruff, mypy)
# 4. Installs package in editable mode
```

### Verify Development Setup

```bash
# Run all quality checks
make check

# Individual commands:
make test      # Run 806 tests with coverage
make lint      # Check code quality (ruff + mypy)
make format    # Auto-format code
```

**Development Commands:**

| Command | Action |
|---------|--------|
| `make install` | Install package for users |
| `make init` | Setup complete dev environment |
| `make test` | Run pytest with coverage (target: 77%+) |
| `make lint` | Check code quality (ruff + mypy) |
| `make format` | Auto-format code with ruff |
| `make check` | Run all quality checks (test + lint) |
| `make clean` | Remove build artifacts |

---

## Read Simulation Setup

External tools required for platform-specific read simulation.

### Illumina Read Simulation

Uses w-Wessim2 pipeline (reseq, bwa, samtools, UCSC tools).

#### Install via Conda

```bash
# Create environment from file
mamba env create -f conda/env_wessim.yaml

# Activate environment
conda activate wessim

# Verify tools
which reseq bwa samtools faToTwoBit pblat
```

#### Manual Installation

Install each tool separately:

1. **reseq** - [https://github.com/schmeing/ReSeq](https://github.com/schmeing/ReSeq)
2. **bwa** - `conda install -c bioconda bwa`
3. **samtools** - `conda install -c bioconda samtools`
4. **UCSC tools** - `conda install -c bioconda ucsc-fatotwobit ucsc-pblat`

Update `config.json` with tool paths:

```json
{
  "tools": {
    "reseq": "/path/to/reseq",
    "bwa": "/path/to/bwa",
    "samtools": "/path/to/samtools",
    "faToTwoBit": "/path/to/faToTwoBit",
    "pblat": "/path/to/pblat"
  }
}
```

---

### Oxford Nanopore Read Simulation

Uses NanoSim with pre-trained error models.

#### Install via Conda

```bash
# Create environment
mamba env create -f conda/env_nanosim.yml

# Activate environment
conda activate nanosim

# Verify installation
which simulator.py minimap2 samtools
```

#### Download Pre-trained Models

NanoSim requires pre-trained error profiles:

```bash
# Example: Human genome R9.4.1 model
wget https://github.com/bcgsc/NanoSim/raw/master/pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training.tar.gz
tar -xzf training.tar.gz

# Update config.json
{
  "nanosim_params": {
    "training_data_path": "/path/to/training"
  }
}
```

---

### PacBio HiFi Read Simulation

Uses pbsim3 for raw reads and CCS for consensus calling.

#### Install via Conda

```bash
# Install pbsim3
conda install -c bioconda pbsim3

# Install PacBio CCS
conda install -c bioconda pbccs

# Install minimap2 and samtools
conda install -c bioconda minimap2 samtools

# Verify installation
which pbsim3 ccs minimap2 samtools
```

---

## Reference Files Setup

Human reference genome and reseq models required for read simulation.

### Automated Download

Use the provided helper script:

```bash
# Download hg38 reference and reseq models
python helpers/install_references.py \
  --assembly hg38 \
  --output-dir reference/

# This downloads:
# - hg38.fa (human genome)
# - hg38.fa.fai (index)
# - reseq error models
```

### Manual Download

#### Human Reference Genome

**hg38 (GRCh38):**
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

**hg19 (GRCh37):**
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
```

Update `config.json`:

```json
{
  "read_simulation": {
    "reference_genome": "/path/to/hg38.fa"
  }
}
```

---

## Shell Completion

Enable autocomplete for faster command-line usage.

### Bash

```bash
_MUCONEUP_COMPLETE=bash_source muconeup > ~/.muconeup-complete.bash
echo 'source ~/.muconeup-complete.bash' >> ~/.bashrc
source ~/.bashrc
```

### Zsh

```bash
_MUCONEUP_COMPLETE=zsh_source muconeup > ~/.muconeup-complete.zsh
echo 'source ~/.muconeup-complete.zsh' >> ~/.zshrc
source ~/.zshrc
```

### Fish

```bash
_MUCONEUP_COMPLETE=fish_source muconeup > ~/.config/fish/completions/muconeup.fish
```

---

## Verify Installation

### Test Core Functionality

```bash
# Check version
muconeup --version

# Validate configuration
muconeup --config config.json simulate --help

# Run quick simulation
muconeup --config config.json simulate \
  --fixed-lengths 20 \
  --out-base test_install
```

### Test Read Simulation

```bash
# Generate test FASTA
muconeup --config config.json simulate --out-base test --fixed-lengths 30

# Test Illumina simulator (if installed)
muconeup --config config.json reads illumina \
  test.001.simulated.fa \
  --coverage 10 \
  --out-base test_reads

# Check output
ls test_reads*
# test_reads_R1.fastq.gz
# test_reads_R2.fastq.gz
# test_reads.illumina.bam
```

---

## Troubleshooting

### Common Issues

**Issue:** `muconeup: command not found`

**Solution:** Ensure installation completed and `~/.local/bin` is in PATH:

```bash
export PATH="$HOME/.local/bin:$PATH"
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
```

---

**Issue:** `ModuleNotFoundError: No module named 'click'`

**Solution:** Reinstall with pip:

```bash
pip install --force-reinstall -e .
```

---

**Issue:** External tool not found (reseq, bwa, etc.)

**Solution:** Install via conda or update `config.json` with correct paths:

```bash
# Check tool location
which reseq

# Update config.json
{
  "tools": {
    "reseq": "/actual/path/to/reseq"
  }
}
```

---

**Issue:** Permission denied on Docker

**Solution:** Add user to docker group:

```bash
sudo usermod -aG docker $USER
# Log out and back in
```

---

## Next Steps

- **[Quick Start Tutorial](quickstart.md)** - Run your first simulation
- **[Core Concepts](concepts.md)** - Understand VNTR simulation fundamentals
- **[Configuration Guide](../reference/configuration.md)** - Customize your setup

---

## Support

- **Issues:** [GitHub Issue Tracker](https://github.com/berntpopp/MucOneUp/issues)
- **Documentation:** [Full Documentation](https://berntpopp.github.io/MucOneUp/)
- **Discussions:** [GitHub Discussions](https://github.com/berntpopp/MucOneUp/discussions)

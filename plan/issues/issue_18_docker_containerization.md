# Issue #18: Docker Containerization (SIMPLIFIED)

## CORRECTION NOTICE
**Original plan proposed 3 separate Dockerfiles + docker-compose, violating KISS principle.**

Simplified approach: **Single Dockerfile, two build targets**

## Problem
Complex conda environment setup (NanoSim Python 3.6, w-Wessim2, external tools). Docker simplifies deployment and ensures reproducible execution environment.

## Simplified Architecture (KISS Principle)

### Single Dockerfile, Two Build Targets:
```
docker/Dockerfile:
  - Target "base": Core app only (haplotype generation)
  - Target "full": Includes all simulators
```

Users choose which image to use:
- **Slim image** (500MB): Haplotype generation only
- **Full image** (2GB): Includes NanoSim, w-Wessim2, all tools

## Implementation

### Single Dockerfile with Multi-Stage Build
```dockerfile
# docker/Dockerfile

# ============================================================================
# Stage 1: Base Dependencies
# ============================================================================
FROM python:3.12-slim AS base-deps

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install uv for fast dependency resolution
RUN pip install --no-cache-dir uv

WORKDIR /app

# ============================================================================
# Stage 2: Slim Image (Core App Only)
# ============================================================================
FROM base-deps AS slim

# Install MucOneUp
COPY pyproject.toml uv.lock ./
RUN uv pip install --system .

# Copy source code
COPY muc_one_up/ /app/muc_one_up/
COPY config.json /app/

# Non-root user
RUN useradd -m -u 1000 muconeup && chown -R muconeup:muconeup /app
USER muconeup

WORKDIR /data
ENTRYPOINT ["muconeup"]
CMD ["--help"]

# ============================================================================
# Stage 3: Full Image (With Simulators)
# ============================================================================
FROM mambaorg/micromamba:1.5.10-bookworm-slim AS full

USER root

# Install system tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    bwa \
    minimap2 \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy MucOneUp installation
COPY --from=slim /usr/local/lib/python3.12 /usr/local/lib/python3.12
COPY --from=slim /usr/local/bin /usr/local/bin
COPY --from=slim /app /app

# Install simulators via conda
USER $MAMBA_USER
COPY conda/env_nanosim.yml /tmp/
COPY conda/env_wessim.yaml /tmp/

RUN micromamba create -y -f /tmp/env_nanosim.yml && \
    micromamba create -y -f /tmp/env_wessim.yaml && \
    micromamba clean --all --yes

# Setup environment activation
ENV PATH="/opt/conda/envs/env_nanosim/bin:/opt/conda/envs/env_wessim/bin:${PATH}"

WORKDIR /data
USER $MAMBA_USER

ENTRYPOINT ["muconeup"]
CMD ["--help"]
```

### Build Script
```bash
#!/bin/bash
# docker/build.sh

set -e

echo "Building MucOneUp Docker images..."

# Build slim image (core app only)
docker build --target slim -t muconeup:slim -f docker/Dockerfile .

# Build full image (with simulators)
docker build --target full -t muconeup:full -f docker/Dockerfile .

# Tag latest as full
docker tag muconeup:full muconeup:latest

echo "✓ Build complete"
echo ""
echo "Images built:"
echo "  muconeup:slim   - Core app (haplotype generation)"
echo "  muconeup:full   - With simulators (Illumina, ONT, PacBio)"
echo "  muconeup:latest - Alias for full"
echo ""
echo "Usage:"
echo "  docker run --rm -v \$(pwd):/data muconeup:slim simulate --help"
echo "  docker run --rm -v \$(pwd):/data muconeup:full reads illumina ..."
```

### .dockerignore
```
.git
.github
__pycache__
*.pyc
.pytest_cache
.mypy_cache
.ruff_cache
htmlcov
.coverage
*.egg-info
output/
tests/
docs/
plan/
.pre-commit-config.yaml
```

## Usage Examples

### Slim Image (Haplotype Generation Only)
```bash
# Build
docker build --target slim -t muconeup:slim .

# Run
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  muconeup:slim \
  --config /app/config.json \
  simulate --out-base /data/test --fixed-lengths 60
```

### Full Image (With Read Simulation)
```bash
# Build
docker build --target full -t muconeup:full .

# Haplotype generation
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  muconeup:full \
  --config /app/config.json \
  simulate --out-base /data/test

# Illumina read simulation
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  muconeup:full \
  --config /app/config.json \
  reads illumina /data/test.fa

# ONT read simulation
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  muconeup:full \
  --config /app/config.json \
  reads ont /data/test.fa
```

### Docker Compose (Optional - For Convenience)
```yaml
# docker-compose.yml (SIMPLE VERSION)
services:
  muconeup:
    image: muconeup:full
    build:
      context: .
      dockerfile: docker/Dockerfile
      target: full
    volumes:
      - ./data:/data
      - ./output:/output
      - ./config.json:/app/config.json:ro
    environment:
      - MUCONEUP_LOG_LEVEL=INFO
```

**Usage**:
```bash
docker-compose run --rm muconeup --config /app/config.json simulate --help
```

## Testing

### Integration Test
```python
def test_docker_slim_image():
    """Test slim image can generate haplotypes."""
    result = subprocess.run([
        "docker", "run", "--rm",
        "-v", f"{tmp_path}:/data",
        "muconeup:slim",
        "--config", "/app/config.json",
        "simulate", "--out-base", "/data/test"
    ], capture_output=True)

    assert result.returncode == 0
    assert (tmp_path / "test.fa").exists()

def test_docker_full_image_illumina():
    """Test full image can simulate Illumina reads."""
    result = subprocess.run([
        "docker", "run", "--rm",
        "-v", f"{tmp_path}:/data",
        "muconeup:full",
        "--config", "/app/config.json",
        "reads", "illumina", "/data/test.fa"
    ], capture_output=True)

    assert result.returncode == 0
    assert (tmp_path / "test.bam").exists()
```

## Documentation

### README.md
```markdown
## Docker Installation

### Quick Start (Slim Image)
```bash
docker pull muconeup:slim
docker run --rm muconeup:slim --help
```

### Full Setup (With Simulators)
```bash
docker pull muconeup:full
docker run --rm -v $(pwd)/data:/data muconeup:full simulate --help
```

### Build from Source
```bash
git clone https://github.com/berntpopp/MucOneUp.git
cd MucOneUp

# Slim image (500MB)
docker build --target slim -t muconeup:slim .

# Full image (2GB)
docker build --target full -t muconeup:full .
```

### Image Comparison
| Image | Size | Includes | Use Case |
|-------|------|----------|----------|
| `muconeup:slim` | ~500MB | Core app only | Haplotype generation |
| `muconeup:full` | ~2GB | + Simulators | Read simulation |
```

## CI/CD Integration

### GitHub Actions
```yaml
# .github/workflows/docker.yml
name: Docker Build

on:
  release:
    types: [published]
  push:
    branches: [main]

jobs:
  docker:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: docker/setup-buildx-action@v3

      - uses: docker/login-action@v3
        with:
          username: ${{ secrets.DOCKERHUB_USERNAME }}
          password: ${{ secrets.DOCKERHUB_TOKEN }}

      - name: Build and push slim
        uses: docker/build-push-action@v6
        with:
          context: .
          file: docker/Dockerfile
          target: slim
          push: true
          tags: muconeup/muconeup:slim

      - name: Build and push full
        uses: docker/build-push-action@v6
        with:
          context: .
          file: docker/Dockerfile
          target: full
          push: true
          tags: |
            muconeup/muconeup:full
            muconeup/muconeup:latest
```

## Files Created

### NEW FILES (3):
- `docker/Dockerfile` (single file with 2 targets)
- `docker/build.sh` (build helper)
- `.dockerignore`
- `docker-compose.yml` (OPTIONAL - simple version)
- `.github/workflows/docker.yml` (CI)

### DO NOT CREATE:
- ❌ `Dockerfile.nanosim` - Unnecessary separate file
- ❌ `Dockerfile.wessim` - Unnecessary separate file
- ❌ Complex docker-compose with profiles - Over-engineering

## Compliance with Programming Principles

✅ **KISS**: Single Dockerfile, simple build targets

✅ **DRY**: Reuses base-deps stage for both images

✅ **Maintainability**: One file to maintain, not three

✅ **User-Friendly**: Clear choice between slim/full images

**Comparison to Original Plan**:
- Original: 3 Dockerfiles + docker-compose + profiles = **Complex**
- Corrected: 1 Dockerfile + 2 targets = **Simple**

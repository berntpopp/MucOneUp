# Issue #18: Docker Containerization - FINAL PLAN (KISS)

**Status:** Implementation Ready
**Priority:** High
**Principle:** **TRUE KISS** - One image, all functionality
**Original Issue:** https://github.com/berntpopp/MucOneUp/issues/18

---

## 🎯 FINAL DESIGN (Ultra-Simplified)

### Single Image Philosophy

**ONE image = muconeup:latest**
- ✅ Core app
- ✅ Illumina simulation (w-Wessim2)
- ✅ ONT simulation (NanoSim)
- ✅ PacBio simulation (pbsim3, pbccs)
- ✅ ALL tools included

**Why?**
- Users don't care about image size (storage is cheap)
- Simpler documentation (one docker pull command)
- No confusion about which image to use
- Single CI/CD workflow
- Single Dockerfile
- Easier maintenance

**Size:** ~2.5GB (acceptable for scientific software with all tools)

---

## 📐 Single Dockerfile (KISS)

```dockerfile
# ============================================================================
# MucOneUp - Single Production Image
# Includes: Core + Illumina + ONT + PacBio simulators
# ============================================================================

FROM mambaorg/micromamba:1.5.10-bookworm-slim

USER root

# System dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        samtools \
        bwa \
        minimap2 \
        git \
        curl \
        ca-certificates \
        python3 \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install MucOneUp Python package
COPY pyproject.toml ./
COPY muc_one_up/ ./muc_one_up/
COPY config.json ./

RUN pip3 install --no-cache-dir --break-system-packages . && \
    pip3 cache purge

# Install ALL conda environments
USER $MAMBA_USER
COPY conda/env_nanosim.yml conda/env_wessim.yaml conda/env_pacbio.yml /tmp/

RUN micromamba create -y -f /tmp/env_nanosim.yml && \
    micromamba create -y -f /tmp/env_wessim.yaml && \
    micromamba create -y -f /tmp/env_pacbio.yml && \
    micromamba clean --all --yes && \
    rm /tmp/*.yml /tmp/*.yaml

# Add all conda env bins to PATH
ENV PATH="/opt/conda/envs/env_nanosim/bin:/opt/conda/envs/env_wessim/bin:/opt/conda/envs/env_pacbio/bin:${PATH}"

# Security
USER $MAMBA_USER
WORKDIR /data

# Metadata
LABEL org.opencontainers.image.source="https://github.com/berntpopp/MucOneUp" \
      org.opencontainers.image.description="MucOneUp VNTR simulator - All-in-one image" \
      org.opencontainers.image.licenses="MIT"

ENTRYPOINT ["muconeup"]
CMD ["--help"]
```

**Simplicity:**
- ✅ No multi-stage complexity
- ✅ No build targets
- ✅ Linear execution flow
- ✅ Easy to understand
- ✅ Easy to maintain

---

## ⚙️ Single GitHub Actions Workflow (KISS)

```yaml
# .github/workflows/docker.yml
name: Docker

on:
  push:
    branches: [main]
  release:
    types: [published]
  pull_request:
    branches: [main]

permissions:
  contents: read
  packages: write
  security-events: write

jobs:
  docker:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      # Extract version
      - name: Get version
        id: version
        run: |
          VERSION=$(grep '__version__' muc_one_up/version.py | cut -d'"' -f2)
          echo "version=$VERSION" >> $GITHUB_OUTPUT

      # Docker setup
      - uses: docker/setup-buildx-action@v3

      - uses: docker/login-action@v3
        if: github.event_name != 'pull_request'
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      # Metadata
      - uses: docker/metadata-action@v5
        id: meta
        with:
          images: ghcr.io/${{ github.repository }}/muconeup
          tags: |
            type=semver,pattern={{version}},value=${{ steps.version.outputs.version }}
            type=semver,pattern={{major}}.{{minor}},value=${{ steps.version.outputs.version }}
            type=edge,branch=main
            type=ref,event=pr
            type=raw,value=latest,enable=${{ github.event_name == 'release' }}

      # Build
      - uses: docker/build-push-action@v6
        with:
          context: .
          file: docker/Dockerfile
          platforms: linux/amd64
          push: ${{ github.event_name != 'pull_request' }}
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha
          cache-to: type=gha,mode=max

      # Security scan
      - uses: aquasecurity/trivy-action@0.30.0
        with:
          image-ref: ghcr.io/${{ github.repository }}/muconeup:${{ steps.version.outputs.version }}
          format: 'sarif'
          output: 'trivy-results.sarif'
          severity: 'CRITICAL,HIGH'
          exit-code: 1

      - uses: github/codeql-action/upload-sarif@v3
        if: always()
        with:
          sarif_file: 'trivy-results.sarif'

      # Test
      - name: Test image
        if: github.event_name != 'pull_request'
        run: |
          docker run --rm ghcr.io/${{ github.repository }}/muconeup:${{ steps.version.outputs.version }} --version
```

**Simplicity:**
- ✅ Single job (no matrix needed)
- ✅ Linear workflow
- ✅ Clear steps
- ✅ No over-engineering

---

## 🔧 Simple docker-compose.yml

```yaml
version: '3.9'

services:
  muconeup:
    image: ghcr.io/berntpopp/muconeup/muconeup:latest
    volumes:
      - ./data:/data
      - ./output:/output
      - ./config.json:/app/config.json:ro
    environment:
      - MUCONEUP_LOG_LEVEL=${LOG_LEVEL:-INFO}
    security_opt:
      - no-new-privileges:true
    cap_drop:
      - ALL
```

**Simplicity:**
- ✅ Single service
- ✅ No YAML anchors needed
- ✅ Copy-paste ready

---

## 📝 Simple Build Script

```bash
#!/usr/bin/env bash
# docker/build.sh

set -euo pipefail

VERSION=$(grep '__version__' muc_one_up/version.py | cut -d'"' -f2)

echo "Building MucOneUp Docker image v${VERSION}"

docker build \
    --tag "muconeup:latest" \
    --tag "muconeup:${VERSION}" \
    --file docker/Dockerfile \
    --progress=plain \
    .

echo "✓ Build complete!"
docker images muconeup

# Test
echo "Testing image..."
docker run --rm muconeup:latest --version
echo "✓ Image working!"
```

**Simplicity:**
- ✅ No loops
- ✅ No arrays
- ✅ Straightforward

---

## 📁 Minimal File Structure

```
MucOneUp/
├── docker/
│   ├── Dockerfile          # Single file, no stages
│   └── build.sh            # Simple build script
├── .dockerignore           # Standard excludes
├── docker-compose.yml      # Single service
└── .github/
    └── workflows/
        └── docker.yml      # Single workflow
```

**4 files total** (vs. 7-10 in over-engineered plans)

---

## 📚 Simple README Section

```markdown
## Docker Installation

### Quick Start

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

### What's Included

One image with everything:
- MucOneUp core
- Illumina simulation (w-Wessim2)
- ONT simulation (NanoSim)
- PacBio simulation (pbsim3, pbccs)

### Docker Compose

```yaml
# docker-compose.yml
version: '3.9'
services:
  muconeup:
    image: ghcr.io/berntpopp/muconeup/muconeup:latest
    volumes:
      - ./data:/data
      - ./config.json:/app/config.json:ro
```

```bash
docker-compose run --rm muconeup --version
```

### Build from Source

```bash
./docker/build.sh
```

Done.
```

**Simplicity:**
- ✅ No image comparison tables
- ✅ No decision paralysis
- ✅ Copy-paste examples
- ✅ Minimal explanation needed

---

## 📊 Implementation Timeline (REVISED)

| Task | Duration |
|------|----------|
| Create Dockerfile | 30 min |
| Create .dockerignore | 10 min |
| Create GitHub Actions workflow | 45 min |
| Create docker-compose.yml | 10 min |
| Create build.sh | 15 min |
| Update README | 30 min |
| Test locally | 30 min |
| **TOTAL** | **3 hours** |

**3 hours vs. 7.5 hours vs. 12-15 hours in previous plans**

---

## ✅ Success Criteria (Simplified)

### Functionality
- [ ] Image builds successfully
- [ ] All tools work (core, Illumina, ONT, PacBio)
- [ ] Version sync works
- [ ] docker-compose works

### Security
- [ ] Trivy scan passes (no HIGH/CRITICAL)
- [ ] Non-root execution
- [ ] Capability dropping in compose

### Quality
- [ ] No duplicated code
- [ ] Clear documentation
- [ ] Easy to maintain
- [ ] No regressions

---

## 🎯 Design Decisions (Rationale)

### Why ONE image?

**Previous thinking (WRONG):**
- "Users might not need all tools"
- "Smaller images are better"
- "Give users choice"

**Correct thinking (KISS):**
- Users don't care about image size (it's 2025, storage is cheap)
- Scientific tools are expected to be large
- Download happens ONCE (cached)
- Simpler mental model = better UX
- Less documentation needed
- Less confusion
- Less testing needed
- Less maintenance

### Why no multi-arch?

**Target audience:** Bioinformaticians using x86-64 servers/workstations
**Reality:** arm64 adds complexity for <1% use case
**KISS:** Build what's needed, not what's possible

### Why no SBOM/provenance in v1?

**v1 goal:** Working Docker image
**v2 goal:** Advanced features
**KISS:** Ship basic version, iterate based on feedback

---

## 🚀 Implementation Order

1. ✅ Create `docker/Dockerfile` (30 min)
2. ✅ Create `.dockerignore` (10 min)
3. ✅ Test local build with `docker/build.sh` (45 min)
4. ✅ Create `.github/workflows/docker.yml` (45 min)
5. ✅ Create `docker-compose.yml` (10 min)
6. ✅ Update README (30 min)
7. ✅ Final testing (30 min)

**Total: 3 hours**

---

## 🔄 Comparison of Plans

| Aspect | Original | Revised | FINAL |
|--------|----------|---------|-------|
| **Images** | 3-4 | 3 | **1** ✅ |
| **Dockerfile stages** | 4 | 4 | **1** ✅ |
| **GitHub Actions jobs** | 5 | 2 | **1** ✅ |
| **Files created** | 10+ | 6 | **4** ✅ |
| **Implementation time** | 12-15h | 7.5h | **3h** ✅ |
| **Maintenance burden** | High | Medium | **Low** ✅ |
| **User confusion** | High | Medium | **None** ✅ |
| **Documentation pages** | 3+ | 2 | **1** ✅ |

---

## 📋 Final Checklist

### Before Implementation
- [ ] Read and understand this plan
- [ ] Verify conda environment files exist
- [ ] Ensure pyproject.toml is correct
- [ ] Check GitHub Actions permissions

### Implementation
- [ ] Create Dockerfile
- [ ] Create .dockerignore
- [ ] Test local build
- [ ] Create GitHub Actions workflow
- [ ] Create docker-compose.yml
- [ ] Update README

### Testing
- [ ] Local build works
- [ ] Image runs successfully
- [ ] All simulators accessible
- [ ] Trivy scan passes
- [ ] docker-compose works
- [ ] CI/CD workflow works

### Deployment
- [ ] Merge to main
- [ ] Verify GHCR publish
- [ ] Test pulling from GHCR
- [ ] Create release (triggers versioned tag)

---

## 💡 Key Principles Applied

### KISS (Keep It Simple, Stupid)
- ✅ **1 image** (not 3, not 4)
- ✅ **1 workflow** (not 3 jobs)
- ✅ **1 Dockerfile** (no multi-stage complexity)
- ✅ **4 total files** (minimal footprint)

### DRY (Don't Repeat Yourself)
- ✅ **Single source of truth** for version
- ✅ **No duplicated code** (no matrix needed with 1 image)
- ✅ **Reusable components** where needed

### YAGNI (You Aren't Gonna Need It)
- ✅ **No SBOM** (not required for v1)
- ✅ **No multi-arch** (not needed)
- ✅ **No multiple images** (users don't need choice)
- ✅ **No weekly scanning** (Dependabot is enough)

### Pragmatism
- ✅ **Ship v1 fast** (3 hours)
- ✅ **Iterate based on feedback**
- ✅ **Add features when requested**
- ✅ **Don't predict future needs**

---

## 🎉 Ready to Implement

This FINAL plan represents the **simplest possible solution** that:
- ✅ Meets all requirements
- ✅ Follows KISS rigorously
- ✅ Eliminates all over-engineering
- ✅ Reduces implementation time to 3 hours
- ✅ Minimizes maintenance burden
- ✅ Maximizes user simplicity

**Next step:** Implement the Dockerfile.

Would you like me to proceed with creating the actual files?

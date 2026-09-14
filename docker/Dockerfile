# syntax=docker/dockerfile:1
# ============================================================================
# MucOneUp - Single Production Image
# Includes: Core + Illumina + ONT + PacBio simulators
# ============================================================================

FROM mambaorg/micromamba:1.5.10-bookworm-slim

USER root

# System dependencies and security upgrades
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        samtools \
        bwa \
        minimap2 \
        git \
        curl \
        ca-certificates \
        python3 \
        python3-pip \
        python3-dev \
        build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Symlink micromamba to mamba so subprocesses calling 'mamba run' work seamlessly
RUN ln -s /usr/bin/micromamba /usr/local/bin/mamba

# Install ALL conda environments FIRST so heavy layers are cached across code changes
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER conda/env_nanosim.yml conda/env_wessim.yaml conda/env_pacbio.yml /tmp/

RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    micromamba create -y -f /tmp/env_nanosim.yml && \
    micromamba create -y -f /tmp/env_wessim.yaml && \
    micromamba create -y -f /tmp/env_pacbio.yml && \
    micromamba clean --all --yes && \
    rm -f /tmp/*.yml /tmp/*.yaml

# Remove unused build and image/font packages from conda environments to harden image and eliminate CVEs
USER root
RUN rm -rf /opt/conda/envs/*/lib/python*/site-packages/{pip,pip-*,setuptools,setuptools-*,wheel,wheel-*,_distutils_hack,distutils-precedence.pth,fonttools,fonttools-*,PIL,Pillow,Pillow-*,pillow,pillow-*}

WORKDIR /app

# Install MucOneUp Python package
COPY pyproject.toml README.md ./
COPY muc_one_up/ ./muc_one_up/
COPY config.json ./

RUN --mount=type=cache,target=/root/.cache/pip \
    pip3 install --no-cache-dir --break-system-packages . biopython && \
    pip3 cache purge

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

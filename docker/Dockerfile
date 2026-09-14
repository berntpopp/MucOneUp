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

RUN --mount=type=cache,target=/home/mambauser/.mamba/pkgs,sharing=locked \
    micromamba create -y -f /tmp/env_nanosim.yml && \
    micromamba create -y -f /tmp/env_wessim.yaml && \
    micromamba create -y -f /tmp/env_pacbio.yml && \
    micromamba install -y -n env_nanosim -c conda-forge --no-deps scikit-learn=1.5.0 && \
    micromamba remove -y -f --no-prune -n env_nanosim pillow fonttools pip wheel setuptools && \
    micromamba remove -y -f --no-prune -n env_wessim pip wheel setuptools && \
    micromamba clean --all --yes && \
    rm -rf /tmp/*.yml /tmp/*.yaml /home/mambauser/.mamba/pkgs/* /opt/conda/pkgs/*

WORKDIR /app

# Install MucOneUp Python package
COPY pyproject.toml README.md ./
COPY muc_one_up/ ./muc_one_up/
COPY config.json ./

USER root
# Install lightweight PIL stubs for matplotlib in env_nanosim (eliminates Pillow CVEs without breaking simulation)
RUN mkdir -p /opt/conda/envs/env_nanosim/lib/python3.9/site-packages/PIL && \
    echo 'class Image: pass' > /opt/conda/envs/env_nanosim/lib/python3.9/site-packages/PIL/Image.py && \
    echo 'class PngInfo: pass' > /opt/conda/envs/env_nanosim/lib/python3.9/site-packages/PIL/PngImagePlugin.py && \
    echo 'from . import Image' > /opt/conda/envs/env_nanosim/lib/python3.9/site-packages/PIL/__init__.py

RUN --mount=type=cache,target=/root/.cache/pip \
    /usr/bin/pip3 install --break-system-packages . biopython && \
    apt-get purge -y python3-pip python3-dev build-essential && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /home/mambauser/.mamba /root/.mamba /tmp/* /var/tmp/* && \
    find /opt/conda -name "__pycache__" -depth -exec rm -rf {} +

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

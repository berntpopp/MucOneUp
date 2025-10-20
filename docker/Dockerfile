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
        python3-pip \
        python3-dev \
        build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install MucOneUp Python package
COPY pyproject.toml README.md ./
COPY muc_one_up/ ./muc_one_up/
COPY config.json ./

RUN pip3 install --no-cache-dir --break-system-packages . biopython && \
    pip3 cache purge

# Install ALL conda environments
USER $MAMBA_USER
COPY conda/env_nanosim.yml conda/env_wessim.yaml conda/env_pacbio.yml /tmp/

RUN micromamba create -y -f /tmp/env_nanosim.yml && \
    micromamba create -y -f /tmp/env_wessim.yaml && \
    micromamba create -y -f /tmp/env_pacbio.yml && \
    micromamba clean --all --yes

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

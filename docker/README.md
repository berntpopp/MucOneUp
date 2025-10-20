# MucOneUp Docker

All-in-one Docker image with MucOneUp + Illumina/ONT/PacBio simulators.

## Installation

```bash
# Pull from GitHub Container Registry
docker pull ghcr.io/berntpopp/muconeup/muconeup:latest

# Verify
docker run --rm ghcr.io/berntpopp/muconeup/muconeup:latest --version
```

## Usage

### Basic Command

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/config.json:/app/config.json:ro \
  ghcr.io/berntpopp/muconeup/muconeup:latest \
  --config /app/config.json \
  simulate --out-base /data/sample
```

### With docker-compose

```bash
# Use docker-compose.yml in repo root
docker-compose run --rm muconeup --config /app/config.json simulate --help
```

### Volume Mounts

- `/data` - Working directory
- `/app/config.json` - Configuration (read-only)

**Example:**
```bash
docker run --rm \
  -v $(pwd)/input:/data \
  -v $(pwd)/output:/output \
  -v $(pwd)/config.json:/app/config.json:ro \
  ghcr.io/berntpopp/muconeup/muconeup:latest \
  --config /app/config.json \
  simulate --out-base /output/result
```

## Build from Source

```bash
# Quick build
./docker/build.sh

# Manual
docker build -t muconeup:latest -f docker/Dockerfile .
```

## Troubleshooting

**Permission errors:** Container runs as UID 1000
```bash
chmod 777 output/  # Or: chown -R 1000:1000 output/
```

**Memory:** Large simulations need more RAM
```bash
docker run --memory=8g --rm ...
```

## Versioning

- `latest` - Stable release
- `0.16.0` - Specific version (recommended for production)
- `edge` - Development (unstable)

**Support:** https://github.com/berntpopp/MucOneUp/issues

#!/usr/bin/env bash
# docker/build.sh - Build MucOneUp Docker image

set -euo pipefail

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Extract version from source
VERSION=$(grep '__version__' muc_one_up/version.py | cut -d'"' -f2)

echo -e "${GREEN}Building MucOneUp Docker image v${VERSION}${NC}"
echo ""

# Build image
echo -e "${YELLOW}Building image...${NC}"
docker build \
    --tag "muconeup:latest" \
    --tag "muconeup:${VERSION}" \
    --file docker/Dockerfile \
    --progress=plain \
    .

echo ""
echo -e "${GREEN}✓ Build complete!${NC}"
docker images muconeup

# Test
echo ""
echo -e "${YELLOW}Testing image...${NC}"
docker run --rm muconeup:latest --version

echo ""
echo -e "${GREEN}✓ Image working!${NC}"
echo ""
echo "Usage:"
echo "  docker run --rm -v \$(pwd)/data:/data muconeup:latest --help"

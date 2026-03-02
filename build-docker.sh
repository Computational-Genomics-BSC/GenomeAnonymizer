#!/bin/bash

# Build script for GenomeAnonymizer Docker image
# This script builds only the Docker image

set -e

# Configuration
IMAGE_NAME="genomeanonymizer"
TAG="latest"
VERSION="1.1.0"

# Check if Docker is available and user has permissions
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not in PATH"
    exit 1
fi

# Test Docker permissions
if ! docker info &> /dev/null; then
    echo "Error: Cannot connect to Docker daemon. This usually means:"
    echo "1. Docker service is not running, or"
    echo "2. Your user doesn't have permission to access Docker"
    echo ""
    echo "To fix this, try one of the following:"
    echo "- Add your user to the docker group: sudo usermod -aG docker \$USER"
    echo "- Then log out and back in, or run: newgrp docker"
    echo "- Or run this script with sudo (not recommended for security)"
    echo ""
    echo "For more info: https://docs.docker.com/engine/install/linux-postinstall/"
    exit 1
fi

echo "Building GenomeAnonymizer Docker image..."

# Build Docker image
docker build -t ${IMAGE_NAME}:${TAG} -t ${IMAGE_NAME}:${VERSION} .

echo "Docker image built successfully!"
echo "Image tags: ${IMAGE_NAME}:${TAG}, ${IMAGE_NAME}:${VERSION}"

# Test the Docker image
echo "Testing Docker image..."
docker run --rm ${IMAGE_NAME}:${TAG} --help

echo "Docker build complete!"
echo ""
echo "Usage examples:"
echo "Basic: docker run --rm ${IMAGE_NAME}:${TAG} --help"
echo "With data: docker run --rm -v /path/to/data:/work ${IMAGE_NAME}:${TAG} [java-options] [app-options]"
echo "Multiple mounts: docker run --rm -v /input:/input -v /output:/output ${IMAGE_NAME}:${TAG} [options]"

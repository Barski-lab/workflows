#!/bin/bash
# deploy_to_hpc.sh
# Script to pull Docker image from DockerHub and convert to Singularity on HPC
# Usage: ./deploy_to_hpc.sh [version]
# If version is not provided, will use 'latest'

set -e  # Exit on any error

# Configuration - modify these variables
IMAGE_NAME="biowardrobe2/scidap-deseq"
SINGULARITY_IMAGE_DIR="/data/barskilab/scidap_server/singularity_images"  # User's actual path
VERSION="${1:-v0.0.30}"  # Use provided version or default to v0.0.30
SINGULARITY_FILENAME="biowardrobe2_scidap-deseq:${VERSION}.sif"  # Match user's naming convention

# Load specific singularity module version
echo "Loading singularity module..."
module load singularity/3.7.0

# Check if Singularity is available
if ! command -v singularity &> /dev/null; then
    echo "Singularity could not be found. Please make sure it's installed or modules are loaded properly."
    exit 1
fi

echo "Using Singularity: $(which singularity)"
echo "Singularity version: $(singularity --version)"

# Create temp directory for download
TEMP_DIR=$(mktemp -d)
cd "${TEMP_DIR}"
echo "Working in temporary directory: ${TEMP_DIR}"

# Remove existing file if it exists
if [ -f "${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}" ]; then
    echo "Removing existing file: ${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
    rm "${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
fi

# Pull the Docker image with Singularity
echo "Pulling ${IMAGE_NAME}:${VERSION} from Docker Hub..."
singularity pull docker://${IMAGE_NAME}:${VERSION}

# Determine pulled filename (singularity usually names it with underscores)
PULLED_FILENAME=$(ls *.sif)
if [ -z "${PULLED_FILENAME}" ]; then
    echo "Failed to pull Singularity image or locate .sif file"
    exit 1
fi

# Ensure destination directory exists
mkdir -p "${SINGULARITY_IMAGE_DIR}"

# Move the image to destination with the desired filename
echo "Moving ${PULLED_FILENAME} to ${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
mv "${PULLED_FILENAME}" "${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"

# Clean up temp directory
cd - > /dev/null
rm -rf "${TEMP_DIR}"

echo "Deployment complete!"
echo "Image available at: ${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"

# Keep a record of the latest deployed version
echo "${VERSION}" > "${SINGULARITY_IMAGE_DIR}/scidap-deseq-version.txt" 
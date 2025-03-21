#!/bin/bash
# check_for_updates.sh
# Script to automatically check for new Docker image versions and pull them with Singularity
# Recommended to set this up as a cron job on the HPC
# Example cron entry (check every day at 2 AM):
# 0 2 * * * /path/to/check_for_updates.sh > /path/to/update_log.txt 2>&1

set -e  # Exit on any error

# Configuration - modify these variables
IMAGE_NAME="biowardrobe2/scidap-deseq"
SINGULARITY_IMAGE_DIR="/data/barskilab/scidap_server/singularity_images"  # User's actual path
VERSION_FILE="${SINGULARITY_IMAGE_DIR}/scidap-deseq-version.txt"

# Log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Load singularity module
log "Loading singularity module..."
module load singularity/3.7.0

# Check if Singularity is available
if ! command -v singularity &> /dev/null; then
    log "Singularity could not be found. Please make sure it's installed or modules are loaded properly."
    exit 1
fi

log "Using Singularity: $(which singularity)"
log "Singularity version: $(singularity --version)"

# Ensure the destination directory exists
mkdir -p "${SINGULARITY_IMAGE_DIR}"

# Get the current version from the version file
if [ -f "${VERSION_FILE}" ]; then
    CURRENT_VERSION=$(cat "${VERSION_FILE}")
    log "Current version: ${CURRENT_VERSION}"
else
    CURRENT_VERSION="none"
    log "No current version found."
fi

# Query DockerHub for the latest tag
log "Checking for latest version on Docker Hub..."
LATEST_VERSION=$(curl -s "https://hub.docker.com/v2/repositories/${IMAGE_NAME}/tags/" | grep -o '"name":"[^"]*' | grep -v latest | head -1 | cut -d'"' -f4)

if [ -z "${LATEST_VERSION}" ]; then
    log "Failed to retrieve latest version information from Docker Hub."
    exit 1
fi

log "Latest version on Docker Hub: ${LATEST_VERSION}"

# Compare versions and pull if different
if [ "${CURRENT_VERSION}" != "${LATEST_VERSION}" ]; then
    log "New version detected. Updating from ${CURRENT_VERSION} to ${LATEST_VERSION}..."
    
    # Create temp directory for download
    TEMP_DIR=$(mktemp -d)
    cd "${TEMP_DIR}"
    log "Working in temporary directory: ${TEMP_DIR}"
    
    # Set the filename based on user's convention
    SINGULARITY_FILENAME="biowardrobe2_scidap-deseq:${LATEST_VERSION}.sif"
    
    # Remove existing file if it exists
    if [ -f "${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}" ]; then
        log "Removing existing file: ${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
        rm "${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
    fi
    
    # Pull the Docker image with Singularity
    log "Pulling ${IMAGE_NAME}:${LATEST_VERSION} from Docker Hub..."
    singularity pull docker://${IMAGE_NAME}:${LATEST_VERSION}
    
    # Determine pulled filename
    PULLED_FILENAME=$(ls *.sif)
    if [ -z "${PULLED_FILENAME}" ]; then
        log "Failed to pull Singularity image or locate .sif file"
        exit 1
    fi
    
    # Move the image to destination with the desired filename
    log "Moving ${PULLED_FILENAME} to ${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
    mv "${PULLED_FILENAME}" "${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
    
    # Clean up temp directory
    cd - > /dev/null
    rm -rf "${TEMP_DIR}"
    
    # Update version file
    echo "${LATEST_VERSION}" > "${VERSION_FILE}"
    
    log "Update complete. Image available at: ${SINGULARITY_IMAGE_DIR}/${SINGULARITY_FILENAME}"
else
    log "Already at the latest version (${CURRENT_VERSION}). No update needed."
fi

log "Check completed successfully." 
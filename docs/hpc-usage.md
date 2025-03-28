# Using Singularity Images on the HPC

This guide explains how to use the Singularity images deployed via the GitHub Actions workflow.

## Available Images

The workflow deploys the following Singularity images to the HPC cluster:

```
/data/barskilab/scidap_server/singularity_images/biowardrobe2_scidap-deseq:v0.0.30.sif
```

The version may change as new versions are deployed. You can check the current version in:
```
/data/barskilab/scidap_server/singularity_images/scidap-deseq-version.txt
```

## Running the Singularity Container

### Interactive Mode

To start an interactive shell session with the container:

```bash
# Load Singularity module
module load singularity/3.7.0

# Start interactive shell
singularity shell /data/barskilab/scidap_server/singularity_images/biowardrobe2_scidap-deseq:v0.0.30.sif
```

### Running Commands

To run a specific command within the container:

```bash
# Load Singularity module
module load singularity/3.7.0

# Run R script
singularity exec /data/barskilab/scidap_server/singularity_images/biowardrobe2_scidap-deseq:v0.0.30.sif Rscript your_analysis.R
```

### Binding Directories

To access directories from the host system within the container:

```bash
# Load Singularity module
module load singularity/3.7.0

# Run with bound directories
singularity exec --bind /path/to/data:/data,/path/to/output:/output \
  /data/barskilab/scidap_server/singularity_images/biowardrobe2_scidap-deseq:v0.0.30.sif \
  Rscript /data/analysis.R
```

## Example Batch Script

Here's an example SLURM batch script to use the container:

```bash
#!/bin/bash
#SBATCH --job-name=deseq_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=deseq_%j.out

# Load Singularity module
module load singularity/3.7.0

# Define input/output paths
INPUT_DIR=/path/to/input
OUTPUT_DIR=/path/to/output
SCRIPT_PATH=/path/to/script.R

# Run analysis in container
singularity exec \
  --bind ${INPUT_DIR}:/input,${OUTPUT_DIR}:/output \
  /data/barskilab/scidap_server/singularity_images/biowardrobe2_scidap-deseq:v0.0.30.sif \
  Rscript ${SCRIPT_PATH}
```

## Troubleshooting

### Common Issues

1. **Module not found**:
   ```
   singularity: command not found
   ```
   Solution: Ensure you've loaded the Singularity module with `module load singularity/3.7.0`

2. **Permission denied**:
   ```
   FATAL: container creation failed: mount /proc/self/fd/3->/var/singularity/mnt/session/rootfs error: permission denied
   ```
   Solution: Check file permissions or use a compute node where you have appropriate permissions

3. **Binding errors**:
   ```
   ERROR: Failed to attach /path/to/data to /data: no such file or directory
   ```
   Solution: Ensure the bind paths exist on both the host and in the container

### Getting Help

For additional help:
1. Check the image documentation in the repository
2. Contact the HPC support team
3. File an issue in the GitHub repository 
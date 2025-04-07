# HPC Deployment Workflow

This repository contains GitHub Actions workflows to automate the deployment of Docker images to HPC clusters using Singularity.

## Workflow Overview

The deployment pipeline consists of two main workflows:

1. **Build and Push Docker Image**
   - Builds a Docker image from the source code
   - Tags the image with the appropriate version
   - Pushes the image to Docker Hub

2. **Deploy to HPC with Singularity**
   - Pulls the Docker image on the HPC compute node
   - Converts it to Singularity format
   - Saves it to a shared location for use by researchers

## How to Use

### Automatic Deployment

The workflow automatically triggers when the "Build and Push Docker Image" workflow completes successfully:

1. The Docker image is built and pushed to Docker Hub
2. The deployment workflow pulls the image to the HPC cluster
3. The Singularity image is created and placed in `/data/barskilab/scidap_server/singularity_images/`

### Manual Deployment

You can also manually trigger the deployment:

1. Go to GitHub Actions → "Deploy to HPC with Singularity" → "Run workflow"
2. Optionally specify a Docker tag (defaults to v0.0.30)
3. Click "Run workflow"

Alternatively, use the GitHub CLI:
```bash
gh workflow run "Deploy to HPC with Singularity" --ref master -f docker_tag=v1.0.0
```

## Self-hosted Runner Setup

This workflow relies on a self-hosted GitHub runner with access to the HPC cluster. The runner must have:

1. SSH access to the HPC cluster (jump host and compute nodes)
2. Proper SSH key configuration (typically `~/.ssh/id_ed25519`)
3. Network connectivity to both GitHub and the HPC cluster

To run the GitHub runner persistently without sudo rights:

```bash
# 1. Create a tmux session that will keep running after you log out
tmux new -s github-runner

# 2. Inside the tmux session, start the runner
cd ~/actions-runner
./run.sh

# 3. Detach from tmux with Ctrl+B, then D
# To reattach later: tmux attach -t github-runner
```

For detailed setup instructions, see [Self-hosted Runner Setup Guide](docs/runner-setup.md).

## Troubleshooting

Common issues and solutions:

- **SSH Connection Failures**: Check SSH keys and network connectivity
- **Singularity Module Issues**: Verify module availability on the compute node
- **Permission Problems**: Ensure file permissions are correct in the destination directory

For detailed logs, check the GitHub Actions output or the `deployment.log` file generated during the run.

## Documentation

- [HPC Usage Guide](docs/hpc-usage.md) - How to use the deployed Singularity images
- [Runner Setup Guide](docs/runner-setup.md) - Detailed instructions for setting up the self-hosted runner

## License

[MIT License](LICENSE)

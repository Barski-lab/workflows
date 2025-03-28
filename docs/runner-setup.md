# Setting Up a Self-hosted GitHub Runner Without Sudo

This guide explains how to set up a self-hosted GitHub Actions runner that persists without requiring sudo privileges, specifically for the HPC deployment workflow.

## Prerequisites

- SSH access to a machine that can reach both GitHub.com and your HPC cluster
- Proper SSH key setup for connecting to the HPC
- `tmux` or `screen` installed (for running the runner persistently)

## Step 1: Download and Configure the Runner

1. In your GitHub repository, go to **Settings** > **Actions** > **Runners** > **New self-hosted runner**

2. Follow the instructions to download and configure the runner:

```bash
# Create a directory
mkdir -p ~/actions-runner && cd ~/actions-runner

# Download the runner package
curl -o actions-runner-linux-x64-2.303.0.tar.gz -L https://github.com/actions/runner/releases/download/v2.303.0/actions-runner-linux-x64-2.303.0.tar.gz

# Extract the installer
tar xzf ./actions-runner-linux-x64-2.303.0.tar.gz

# Configure the runner
./config.sh --url https://github.com/YOUR_USERNAME/YOUR_REPO --token YOUR_TOKEN
```

## Step 2: Set Up a Persistent Session

Since we can't use systemd without sudo rights, we'll use tmux to keep the runner running:

```bash
# Install tmux if not already available (you may need to use a package manager available to you)
# If you can't install it, check if it's already available:
which tmux

# Create a new tmux session
tmux new -s github-runner

# Inside the tmux session, start the runner
cd ~/actions-runner
./run.sh

# Detach from the tmux session with:
# Press Ctrl+B, then press D
```

To reattach to the session later:
```bash
tmux attach -t github-runner
```

## Step 3: Verify SSH Configuration

Ensure your SSH configuration is correctly set up:

```bash
# Check if your SSH key exists
ls -la ~/.ssh/id_ed25519

# If you need to generate a key:
ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519

# Test connection to the jump host
ssh bmicluster "hostname && whoami"

# Test connection to the compute node
ssh bmicluster-compute "hostname && whoami"
```

## Step 4: Test the Runner

1. Push a small change to your repository or manually trigger a workflow
2. Go to the Actions tab in your repository to see if the runner picks up the job
3. Check the job logs for any errors

## Maintenance

### Restarting After Reboot

If the machine reboots, you'll need to restart the runner:

```bash
# Reattach to existing session (it might be gone after reboot)
tmux attach -t github-runner

# If the session doesn't exist, create a new one
tmux new -s github-runner

# Start the runner
cd ~/actions-runner
./run.sh
```

### Updating the Runner

To update the runner software:

```bash
# Stop any running runner by pressing Ctrl+C in the tmux session
# Then run:
cd ~/actions-runner
./config.sh --url https://github.com/YOUR_USERNAME/YOUR_REPO --token YOUR_TOKEN
```

## Troubleshooting

### Runner Not Connecting

If the runner can't connect to GitHub:

1. Check internet connectivity
2. Verify firewall settings allow outbound HTTPS connections
3. Check GitHub status page

### Jobs Failing with SSH Errors

If jobs fail with SSH connection errors:

1. Verify SSH keys are set up correctly
2. Check if SSH config has the correct hostnames and usernames
3. Try connecting manually to verify credentials work 
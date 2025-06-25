#!/bin/bash

# FilTar Development Environment Setup Script
# This script sets up the conda environment and installs dependencies

set -e

echo "Setting up FilTar development environment..."

# Initialize conda for the current user
conda init bash
source ~/.bashrc

# Configure conda to not auto-activate base environment
conda config --set auto_activate_base false

# Update conda
conda update -y -n base -c defaults conda

# Install system dependencies that are commonly needed
sudo apt-get update
sudo apt-get install -y gzip wget curl

sudo ln -s /bin/tar /bin/gtar || true

# Create the conda environment from environment.yml
echo "Creating conda environment from environment.yml..."
conda env create -f environment.yml

# Reinitialize conda to ensure conda activate function is available
eval "$(conda shell.bash hook)"

# Activate the environment
echo "Activating test-environment..."
conda activate test-environment

# Verify installations
echo "Verifying installations..."
python --version
conda --version
snakemake --version
R --version

# Create workspace settings for Python interpreter
echo "Configuring workspace settings..."
mkdir -p .vscode
cat > .vscode/settings.json << EOF
{
    "python.defaultInterpreterPath": "/opt/conda/envs/test-environment/bin/python",
    "python.terminal.activateEnvironment": true
}
EOF

echo "Development environment setup complete!"
echo "To activate the environment, run: conda activate test-environment"

# Configure shell to auto-activate test-environment for new terminals
echo "" >> ~/.bashrc
echo "# Auto-activate test-environment for FilTar development" >> ~/.bashrc
echo "if [[ \$CONDA_DEFAULT_ENV != 'test-environment' ]]; then" >> ~/.bashrc
echo "    conda activate test-environment 2>/dev/null || true" >> ~/.bashrc
echo "fi" >> ~/.bashrc
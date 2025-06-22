#!/bin/bash

# FilTar Development Environment Setup Script
# This script sets up the conda environment and installs dependencies

set -e

echo "Setting up FilTar development environment..."

# Initialize conda for the current user
conda init bash
source ~/.bashrc

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

# Install R dependencies
echo "Installing R dependencies..."
Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("TBradley27/filtar_R")'
Rscript -e 'install.packages("testthat", repos="https://cloud.r-project.org")'

# Verify installations
echo "Verifying installations..."
python --version
conda --version
snakemake --version
R --version

echo "Development environment setup complete!"
echo "To activate the environment, run: conda activate test-environment"
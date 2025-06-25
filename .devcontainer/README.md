# FilTar Development Container

This directory contains the development container configuration for the FilTar project. The devcontainer provides a consistent development environment that can be easily spun up in GitHub Codespaces, VS Code Remote Containers, or other devcontainer-compatible tools.

## Features

### Essential Requirements ✓
- **Conda Environment**: Pre-configured with conda for package management
- **Snakemake 8**: Snakemake 8.30.0 installed for workflow management

### Desired Requirements ✓
- **R Support**: R development environment with syntax highlighting and debugging extensions
- **Snakemake Extensions**: Syntax highlighting and language support for Snakemake workflows
- **GitHub Integration**: GitHub Actions and Pull Request extensions for seamless GitHub workflow
- **Python Support**: Full Python development environment with linting and debugging capabilities

## Usage

### GitHub Codespaces
1. Open the repository in GitHub
2. Click the "Code" button and select "Codespaces"
3. Click "Create codespace on main" (or your preferred branch)
4. Wait for the container to build and the environment to be set up

### VS Code Remote Containers
1. Install the "Remote - Containers" extension in VS Code
2. Open the repository in VS Code
3. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac) and select "Remote-Containers: Reopen in Container"
4. Wait for the container to build and the environment to be set up

## Environment Setup

The container automatically:
1. Uses the Rocker R 4.2 base image with R pre-installed
2. Installs essential R packages (remotes, testthat, languageserver) via devcontainer features
3. Creates a conda environment named `test-environment` from `environment.yml`
4. Installs Snakemake 8.30.0 and Python 3.12
5. Installs the FilTar R package from GitHub
6. Configures VS Code with appropriate extensions for R, Python, Snakemake, and GitHub workflows

## Activating the Environment

The conda environment `test-environment` is automatically activated when you open a new terminal. If you need to manually activate it, use:
```bash
conda activate test-environment
```

## Included Extensions

- **Python**: `ms-python.python`, `ms-python.pylint`, `ms-python.flake8`, `ms-python.autopep8`
- **R**: `reditorsupport.r`, `rdebugger.r-debugger`, `posit.air-vscode`
- **Snakemake**: `snakemake-lang.snakemake-lang`, `tfehlmann.snakefmt`
- **GitHub**: `github.vscode-github-actions`, `github.vscode-pull-request-github`
- **Additional**: JSON, YAML support

## Troubleshooting

### Common Issues and Solutions

**R Terminal Issues**: If R terminals cannot attach, the devcontainer automatically configures `r.rterm.linux` and `r.rpath.linux` to point to the R installation provided by the Rocker base image.

**Python Interpreter Issues**: The devcontainer creates workspace-specific settings to ensure the correct Python interpreter from the conda environment is used. If VS Code doesn't detect it automatically, check the bottom-left status bar and select the correct interpreter.

**GitHub Extension Issues**: If GitHub-related extensions show commit errors during initial setup, these typically resolve once the environment is fully built and VS Code restarts.

### Manual Troubleshooting Steps
If you encounter issues with the environment setup:
1. Check the terminal output during container creation
2. Manually run the setup script: `bash .devcontainer/setup.sh`
3. Verify the conda environment: `conda env list`
4. Check Snakemake installation: `snakemake --version`
5. Verify R installation: `R --version`
6. Check Python interpreter: look for `/opt/conda/envs/test-environment/bin/python` in VS Code status bar
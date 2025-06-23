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
1. Creates a conda environment named `test-environment` from `environment.yml`
2. Installs Snakemake 8.30.0 and Python 3.12
3. Installs R dependencies from `r-requirements.txt`
4. Configures VS Code with appropriate extensions for R, Python, Snakemake, and GitHub workflows

## Activating the Environment

Once the container is running, you can activate the conda environment with:
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

If you encounter issues with the environment setup:
1. Check the terminal output during container creation
2. Manually run the setup script: `bash .devcontainer/setup.sh`
3. Verify the conda environment: `conda env list`
4. Check Snakemake installation: `snakemake --version`
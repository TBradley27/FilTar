# Contributing: CI/CD Setup for FilTar

This document describes the Continuous Integration and Continuous Deployment (CI/CD) setup for the FilTar repository. It is intended for contributors and maintainers who want to understand or modify the automated workflows that ensure code quality and reproducibility.

## Overview

FilTar uses GitHub Actions for CI/CD. The main workflow is defined in `.github/workflows/ci.yml` and is triggered on pushes and pull requests to the `master` and `use_github_actions` branches.

## Key Features

- **Python & Snakemake Environment:**
  - Managed with Conda using an `environment.yml` file.
  - The environment is cached to speed up builds.
  - Environment is created or updated automatically in the workflow.

- **R Environment:**
  - R is set up using the official `r-lib/actions/setup-r` GitHub Action.
  - R package dependencies are listed in `r-requirements.txt` and installed via Rscript commands in the workflow.
  - R package library is cached using GitHub Actions cache to avoid repeated downloads and installations.

- **System Dependencies:**
  - Required system packages (e.g., gzip, gtar) are installed at the start of the workflow.

- **Testing:**
  - The workflow runs Snakemake pipelines and R tests to ensure code correctness.

## Workflow Steps

1. **Checkout repository**
2. **Set up R and cache R packages**
3. **Set up Miniconda and cache Conda environment**
4. **Install system dependencies**
5. **Create or update Conda environment**
6. **Install R dependencies**
7. **Run Snakemake pipelines and R tests**

## Caching

- **Conda cache:**
  - Caches Conda packages and environments based on the hash of `environment.yml`.
  - If the environment file changes, the cache is invalidated and rebuilt.
- **R package cache:**
  - Caches the R library directory based on the hash of `r-requirements.txt`.
  - If R requirements change, the cache is invalidated and rebuilt.

## Modifying the CI/CD Setup

- To add or update Python/Snakemake dependencies, edit `environment.yml`.
- To add or update R dependencies, edit `r-requirements.txt` and the relevant install commands in `.github/workflows/ci.yml`.
- To add system dependencies, update the relevant step in `.github/workflows/ci.yml`.
- For workflow logic changes, edit `.github/workflows/ci.yml` directly.

## Further Reading

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Conda Documentation](https://docs.conda.io/en/latest/)
- [Snakemake Documentation](https://snakemake.readthedocs.io)
- [R-lib/actions/setup-r](https://github.com/r-lib/actions/tree/v2/setup-r)

---

For any questions or suggestions, please open an issue or pull request.

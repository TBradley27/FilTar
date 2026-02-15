[![CI](https://github.com/TBradley27/FilTar/actions/workflows/ci.yml/badge.svg)](https://github.com/TBradley27/FilTar/actions/workflows/ci.yml)
[![GitHub release](https://img.shields.io/github/release/TBradley27/FilTar.svg)](https://GitHub.com/TBradley27/FilTar/releases/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# FilTar

FilTar is a tool to integrate RNA-Seq data to pre-existing miRNA target prediction workflows in order to increase prediction accuracy.

It achieves this by:

1. Removing transcripts which are not expressed or poorly expressed for a given cell type or tissue
2. Generating 3'UTR annotations specific to a given cell type or tissue

It also operates as a fully functional wrapper around the pre-existing TargetScan7 and miRanda target prediction workflows.

## Installation

Instructions on how to install FilTar can be found at the following location: https://tbradley27.github.io/FilTar/

## Basic Usage

FilTar can be used by following 2 steps:

1. Specify the options you would like to use to run FilTar by editing `config/basic.yaml`.
2. Run the following command:
```
snakemake --use-conda --cores $N target_predictions.txt
```

After running the command, all target predictions are contained inside `target_predictions.txt`.

The following video presents a concise demonstration of basic FilTar usage:

https://www.youtube.com/watch?v=Xhl-nsg7_xo

More detailed instructions can be found inside the full documentation: https://tbradley27.github.io/FilTar/

## Repository Structure

FilTar follows Snakemake best practices for workflow deployment:

```
FilTar/
├── workflow/              # Main workflow directory
│   ├── Snakefile         # Core workflow definition
│   ├── rules/            # Workflow modules (data download, analysis, etc.)
│   ├── scripts/          # All workflow scripts (Python, R, shell)
│   └── envs/             # Conda environment definitions
├── config/               # Configuration files
├── data/                 # Input data (downloaded/user-provided)
└── results/              # Analysis outputs
```

The workflow structure ensures:
- **Modularity**: Clear separation of different analysis components
- **Reproducibility**: Consolidated environment and script management
- **Maintainability**: Organized codebase following established conventions
- **Standards Compliance**: Follows Snakemake's automatic workflow detection

## Containerisation

FilTar provides a development container configuration that offers a consistent, pre-configured development environment for seamless collaboration and development. The containerised environment includes all necessary dependencies (R, Python, Snakemake, conda) and VS Code extensions optimised for FilTar development.

This approach is particularly useful for:
- **GitHub Codespaces**: One-click cloud development environment
- **VS Code Remote Containers**: Local containerised development
- **Consistent environments**: Eliminates "works on my machine" issues

To get started with the containerised environment, simply open the repository in GitHub Codespaces or use VS Code's Remote Containers extension. For detailed setup instructions, troubleshooting, and feature documentation, see the [development container README](.devcontainer/README.md).

## CI/CD and Contributing

For details on the automated CI/CD setup, see the [contributing guide](docs/contributing.md).

This document explains how FilTar uses GitHub Actions for continuous integration and deployment, including environment management, caching, and automated testing. Contributors should review it before making changes to the workflow or dependencies.

## Publication

The article describing FilTar can be found in Volume 36, Issue 8 (pages 2410-2416) of the *Bioinformatics* journal published by Oxford University Press. An online, open access version of the article is available [here](https://doi.org/10.1093/bioinformatics/btaa007 "FilTar Bioinformatics article").

__DOI:__ [10.1093/bioinformatics/btaa007](https://doi.org/10.1093/bioinformatics/btaa007)

__PMCID:__ [PMC7178423](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7178423/)

__PMID:__ 31930382 

The default method of citing the article is to use the following:

> Thomas Bradley, Simon Moxon, FilTar: using RNA-Seq data to improve microRNA target prediction accuracy in animals, Bioinformatics, Volume 36, Issue 8, 15 April 2020, Pages 2410–2416, https://doi.org/10.1093/bioinformatics/btaa007

## Getting Help
If you would like help using FilTar, and are having issues not related to a bug, please raise a new item in the 'discussions' section of this repo

## Reporting bugs, suggested enhancements or any other issues

The issues page of this repository is the best place to post this.

## Contributions and Acknowledgements

Simox Moxon came up with the original idea and project proposal for FilTar. The FilTar concept was extended and developed further between Simon Moxon and Thomas Bradley through the course of the latter's BBSRC (Biotechnology and Biological Sciences Research Council) PhD Studentship on the Norwich Research Park (NRP) Bioscience Doctoral Training Partnership (DTP) programme, when Thomas Bradley worked under the primary supervision of Simon Moxon initially predominantly at the Earlham Institute, and then later predominantly at the School of Biological Sciences, University of East Anglia.

Full acknowledgements can be found within the acknowledgements section of the FilTar publication.

# Workflow Structure

This directory contains the FilTar workflow following Snakemake best practices for deployment.

## Directory Structure

```
workflow/
├── Snakefile              # Main workflow entry point
├── rules/                 # Individual workflow modules (formerly modules/)
│   ├── data_download/     # Data download rules
│   ├── trim_reads/        # Read trimming rules
│   ├── quant_reads/       # Read quantification rules
│   ├── target_prediction/ # miRNA target prediction rules
│   └── ...                # Other workflow modules
├── scripts/               # All workflow scripts (Python, R, shell)
├── envs/                  # All conda environment files
└── notebooks/             # Jupyter notebooks (if any)
```

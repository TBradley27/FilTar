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

## Key Changes

1. **Consolidated Environment Files**: All conda environment files have been moved from scattered `envs/` directories within modules to a single `workflow/envs/` directory.

2. **Centralized Scripts**: All scripts (Python, R, shell) have been moved to `workflow/scripts/` for better organization and reusability.

3. **Modules → Rules**: The `modules/` directory has been renamed to `rules/` following Snakemake conventions.

4. **Automatic Detection**: Snakemake automatically detects and uses this workflow directory, so existing usage patterns remain unchanged:
   ```bash
   snakemake --use-conda --cores $N target_predictions.txt
   ```

## Environment Files

Environment file conflicts have been resolved:
- `bedtools-general.yaml`: Used by with_reannotation module
- `bedtools-v2-27-1.yaml`: Used by get_utr_and_cds and miRanda modules
- All other environment files are consolidated and deduplicated

## Path Resolution

All paths in the workflow are resolved relative to the repository root, ensuring consistent access to data, results, and configuration files.
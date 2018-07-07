#!/bin/bash

#BSUB -n 16              # number of processors required
#BSUB -J snakemake               # name of job
#BSUB -oo snakemake.out          # log file for standard output (overwrites)
#BSUB -eo snakemake.err          # log file for standard error (overwrites)
#BSUB -q moxon
#BSUB -M 4096           # upper memory limit in Mb

snakemake --cores 16 --use-conda $1

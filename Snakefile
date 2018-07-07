#!/bin/bash

configfile: "config.yaml"

include: "sub_snakemake/data_download/ENA/Snakefile"
include: "sub_snakemake/qc/Snakefile"
include: "sub_snakemake/trim_reads/trim_galore/Snakefile"
include: "sub_snakemake/map_reads/Snakefile"
include: "sub_snakemake/APAtrap/Snakefile"
include: "sub_snakemake/reannotate_3utrs/Snakefile"

rule all:
     input: expand("results/{accession}_trimmed.fq.gz", accession=config['all_single_end'])

# expand("data/{accession}.fastq.gz", accession=config['projects']['mouse']['single-end']['PRJNA143627']) #expand("results/{accession}_trimmed.fq.gz", accession=config['all_single_end'])


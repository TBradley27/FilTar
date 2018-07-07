#!/bin/bash

include: "sub_snakemake/data_download/Snakefile"
include: "sub_snakemake/qc/Snakefile"
include: "sub_snakemake/trim_reads/Snakefile"
include: "sub_snakemake/map_reads/Snakefile"
include: "sub_snakemake/APAtrap/Snakefile"
include: "sub_snakemake/reannotate_3utrs/Snakefile"

configfile: "config.yaml"

rule all:
     input: expand("results/{accession}.bam", accession=config['all_brain_runs'])


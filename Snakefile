#!/bin/bash

include: "sub_snakemake/data_download/Snakefile", "sub_snakemake/qc/Snakefile", "sub_snakemake/trim_reads/Snakefile",
         "sub_snakemake/map_reads/Snakefile", "sub_snakemake/APAtrap/Snakefile", "sub_snakemake/reannotate_3utrs/Snakefile"

configfile: "config.yaml"

rule all:
     input: expand("results/{accession}.bam", accession=config['all_brain_runs'])


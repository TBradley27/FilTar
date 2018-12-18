configfile: "config/config.yaml"
configfile: "config/species_basic.yaml"

include: "sub_snakemake/data_download/ENA/Snakefile"
include: "sub_snakemake/data_download/Snakefile"
include: "sub_snakemake/quant_reads/Snakefile"
include: "sub_snakemake/trim_reads/trim_galore/Snakefile"
include: "sub_snakemake/map_reads/hisat2/Snakefile"
include: "sub_snakemake/APAtrap/Snakefile"
include: "sub_snakemake/reannotate_3utrs/Snakefile"
include: "sub_snakemake/quant_reads/salmon/Snakefile"
include: "sub_snakemake/mirna/Snakefile"
include: "sub_snakemake/target_prediction/targetscan/Snakefile"
include: "sub_snakemake/get_utr_and_cds/Snakefile"
include: "sub_snakemake/canonical_targets/Snakefile"

rule all:
     input: expand("results/targets/canonical/{species}_chr{chrom}_msa.contextpp.tsv", species='hsa', chrom='Y')

configfile: "config/accession_mappings.yaml"
configfile: "config/matedness.yaml"
configfile: "config/species_basic.yaml"
configfile: "config/species_sequencing.yaml"

include: "sub_snakemake/data_download/ENA/Snakefile"
#include: "sub_snakemake/qc/Snakefile"
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
include: "sub_snakemake/target_prediction/miRanda/Snakefile"
include: "sub_snakemake/create_tables/SQLite/Snakefile"
include: "sub_snakemake/upload_to_tables/SQLite/Snakefile"
include: "sub_snakemake/upload_to_tables/Snakefile"
include: "sub_snakemake/profiling/Snakefile"
include: "sub_snakemake/no_reannotation/Snakefile"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z-9]{1,2}"

rule all:
     input: expand("results/targets/canonical/{species}_chr{chrom}_msa.contextpp.tsv", species='hsa', chrom='Y')

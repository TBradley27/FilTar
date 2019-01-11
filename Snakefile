configfile: "config/config.yaml"
configfile: "config/species_basic.yaml"

include: "sub_snakemake/without_reannotation/Snakefile"
include: "sub_snakemake/data_download/ENA/Snakefile"
include: "sub_snakemake/data_download/Snakefile"
include: "sub_snakemake/trim_reads/trim_galore/Snakefile"
#include: "sub_snakemake/with_reannotation/map_reads/hisat2/Snakefile"
include: "sub_snakemake/quant_reads/salmon/Snakefile"
include: "sub_snakemake/mirna/Snakefile"
include: "sub_snakemake/target_prediction/targetscan/Snakefile"
include: "sub_snakemake/get_utr_and_cds/without_conservation/Snakefile" # for conservation information substitute 'no_conservation' for 'with_conservation'
include: "sub_snakemake/target_prediction/miRanda/Snakefile"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,5}",
    feature="(3UTR|CDS)",
    ensembl_release="[0-9]{2,3}",
    genus_species="[A-Z][a-z]+_[a-z]+"

rule all:
     input: #"results/targets/mmu/liver.contextpp.tsv" 

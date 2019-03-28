import re

configfile: "config/config.yaml"
configfile: "config/species_basic.yaml"
configfile: "config/validation.yaml"

include: "sub_snakemake/with_reannotation/Snakefile"
include: "sub_snakemake/with_reannotation/map_reads/hisat2/Snakefile"
include: "sub_snakemake/without_reannotation/Snakefile"

include: "sub_snakemake/target_prediction/miRanda/Snakefile"

if config['sequence_data_source'] == 'ENA':
	include: "sub_snakemake/data_download/ENA/Snakefile"
elif config['sequence_data_source'] == 'SRA':
	include: "sub_snakemake/data_download/SRAtoolkit/Snakefile"
elif config['sequence_data_source'] == 'N/A':
	pass
else:
	raise Exception("\nPlease enter a value of either 'ENA' or 'SRA' or 'N/A' for the 'sequence_data_source' key. Default values can be set in config/config.yaml\n")

for transcript in list(config['transcripts']):
	if re.match('^ENS[A-Z]+[0-9]+.[1-9]$',transcript):
		pass
	else:
		raise Exception('\nInvalid transcript identifier "{}". Identifiers must adhere to official Ensembl identifier patterns e.g. "ENSMUST00000189888.6". Please revise.\n'.format(transcript))

include: "sub_snakemake/qc/Snakefile"
include: "sub_snakemake/data_download/Snakefile"
include: "sub_snakemake/trim_reads/trim_galore/Snakefile"
include: "sub_snakemake/quant_reads/kallisto/Snakefile"
include: "sub_snakemake/quant_reads/salmon/Snakefile"
include: "sub_snakemake/mirna/Snakefile"
include: "sub_snakemake/target_prediction/targetscan/Snakefile"
include: "sub_snakemake/cumulative_plots/Snakefile"
include: "sub_snakemake/view_bams/Snakefile"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,5}",
    feature="(3UTR|CDS)",
    ensembl_release="[0-9]{2,3}",
    genus_species="[A-Z][a-z]+_[a-z]+"

rule all:
     input:
             "results/plots/supplementary_data3.pdf","results/plots/supplementary_data4.pdf",
             "results/plots/supplementary_data5.pdf" 

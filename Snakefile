import re

configfile: "config/config.yaml"
configfile: "config/species_basic.yaml"
configfile: "config/species_sequencing.yaml"

if config['reannotation'] == True:
	include: "sub_snakemake/with_reannotation/Snakefile"
	include: "sub_snakemake/with_reannotation/map_reads/hisat2/Snakefile"
elif config['reannotation'] == False:
	include: "sub_snakemake/without_reannotation/Snakefile"
else:
	raise Exception("\nPlease enter a value of either 'True' or 'False' for the 'reannotation' key. Default values can be set in config/config.yaml\n")

if config['conservation'] == True:
	include: "sub_snakemake/get_utr_and_cds/with_conservation/Snakefile"
elif config['conservation'] == False:
	include: "sub_snakemake/get_utr_and_cds/without_conservation/Snakefile"
	include: "sub_snakemake/target_prediction/miRanda/Snakefile"
else:
	raise Exception("\nPlease enter a value of either 'True' or 'False' for the 'conservation' key. Default values can be set in config/config.yaml\n")

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

include: "sub_snakemake/data_download/Snakefile"
include: "sub_snakemake/trim_reads/trim_galore/Snakefile"
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
include: "sub_snakemake/get_utr_and_cds/no_conservation/Snakefile" # for conservation information substitute 'no_conservation' for 'with_conservation'
include: "sub_snakemake/get_utr_and_cds/without_conservation/Snakefile" # for conservation information substitute 'no_conservation' for 'with_conservation'
include: "sub_snakemake/target_prediction/miRanda/Snakefile"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,5}",
    feature="(3UTR|CDS)",
    ensembl_release="[0-9]{2,3}",
    genus_species="[A-Z][a-z]+_[a-z]+"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,2}",
    feature="(3UTR|CDS)"

rule all:
     input: expand("results/targets/canonical/{species}_chr{chrom}_msa.contextpp.tsv", species='hsa', chrom='Y')

#FilTar: Integrating RNA-Seq data to improve microRNA target prediction accuracy in animals
#Copyright (C) 2019 Thomas Bradley
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

import re

configfile: "config/basic.yaml"
configfile: "config/species.yaml"
configfile: "config/accession_mappings.yaml"

if config['reannotation'] == True:
	include: "modules/with_reannotation/Snakefile"
	include: "modules/with_reannotation/map_reads/hisat2/Snakefile"
elif config['reannotation'] == False:
	include: "modules/without_reannotation/Snakefile"
else:
	raise Exception("\nPlease enter a value of either 'True' or 'False' for the 'reannotation' key. Default values can be set in config/config.yaml\n")

if config['conservation'] == True:
	include: "modules/get_utr_and_cds/with_conservation/Snakefile"
elif config['conservation'] == False:
	include: "modules/get_utr_and_cds/without_conservation/Snakefile"
	include: "modules/target_prediction/miRanda/Snakefile"
else:
	raise Exception("\nPlease enter a value of either 'True' or 'False' for the 'conservation' key. Default values can be set in config/config.yaml\n")

if config['sequence_data_source'] == 'ENA':
	include: "modules/data_download/ENA/Snakefile"
elif config['sequence_data_source'] == 'SRA':
	include: "modules/data_download/SRAtoolkit/Snakefile"
elif config['sequence_data_source'] == 'N/A':
	pass
else:
	raise Exception("\nPlease enter a value of either 'ENA' or 'SRA' or 'N/A' for the 'sequence_data_source' key. Default values can be set in config/config.yaml\n")

for transcript in list(config['transcripts']):
	if re.match('^ENS[A-Z]+[0-9]+.[1-9]{1,2}$',transcript):
		pass
	else:
		raise Exception('\nInvalid transcript identifier "{}". Identifiers must adhere to official Ensembl identifier patterns e.g. "ENSMUST00000189888.6". Please revise.\n'.format(transcript))

include: "modules/data_download/Snakefile"
include: "modules/trim_reads/trim_galore/Snakefile"
include: "modules/quant_reads/salmon/Snakefile"
include: "modules/mirna/Snakefile"
include: "modules/target_prediction/targetscan/Snakefile"

include: "modules/create_tables/SQLite/Snakefile"
include: "modules/upload_to_tables/Snakefile"
include: "modules/upload_to_tables/SQLite/Snakefile"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,5}",
    feature="(3UTR|CDS)",
    ensembl_release="[0-9]{2,3}",
    genus_species="[A-Z][a-z]+_[a-z]+"

rule all:
	input: "results/gene_table_loaded_hsa.txt","results/gene_table_loaded_mmu.txt",
               "results/gene_species_table_loaded_hsa.txt","results/gene_species_table_loaded_mmu.txt",
               "results/mRNA_table_loaded_hsa.txt","results/mRNA_table_loaded_mmu.txt",
               "results/Tissues_table_loaded_hsa.txt","results/Tissues_table_loaded_mmu.txt",
               "results/samples_table_loaded_hsa.txt","results/samples_table_loaded_mmu.txt",
               "results/runs_table_loaded_hsa.txt","results/runs_table_loaded_mmu.txt",
               "results/mirnas_loaded.txt"

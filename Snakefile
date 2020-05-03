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
configfile: "config/validation.yaml"

include: "modules/with_reannotation/Snakefile"
include: "modules/with_reannotation/map_reads/hisat2/Snakefile"
include: "modules/without_reannotation/Snakefile"

include: "modules/target_prediction/miRanda/Snakefile"

if config['sequence_data_source'] == 'ENA':
	include: "modules/data_download/ENA/Snakefile"
elif config['sequence_data_source'] == 'SRA':
	include: "modules/data_download/SRAtoolkit/Snakefile"
elif config['sequence_data_source'] == 'N/A':
	pass
else:
	raise Exception("\nPlease enter a value of either 'ENA' or 'SRA' or 'N/A' for the 'sequence_data_source' key. Default values can be set in config/config.yaml\n")

for transcript in list(config['transcripts']):
	if re.match('^ENS[A-Z]+[0-9]+.[1-9]$',transcript):
		pass
	else:
		raise Exception('\nInvalid transcript identifier "{}". Identifiers must adhere to official Ensembl identifier patterns e.g. "ENSMUST00000189888.6". Please revise.\n'.format(transcript))

include: "modules/qc/Snakefile"
include: "modules/data_download/Snakefile"
include: "modules/trim_reads/trim_galore/Snakefile"
include: "modules/quant_reads/kallisto/Snakefile"
include: "modules/quant_reads/salmon/Snakefile"
include: "modules/mirna/Snakefile"
include: "modules/target_prediction/targetscan/Snakefile"
include: "modules/cumulative_plots/Snakefile"
include: "modules/view_bams/Snakefile"

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,5}",
    feature="(3UTR|CDS)",
    ensembl_release="[0-9]{2,3}",
    genus_species="[A-Z][a-z]+_[a-z]+"

rule all:
        input:
                #"results/plots/supplementary_data1.pdf","results/plots/supplementary_data3.pdf",
                #"results/plots/supplementary_data4.pdf","results/plots/supplementary_data4.pdf",
                #'results/plots/figure1.png','results/plots/supplementary_figure1.png',
                #'results/plots/figure2.png','results/plots/figure3.png',
                #'results/plots/supplementary_figure7.png','results/plots/supplementary_figure2.png',
                #"results/plots/supplementary_figure3a.png","results/plots/supplementary_figure3b.png",
                #'results/plots/supplementary_table4.tsv','results/plots/figure4.png',
                #'results/plots/supplementary_figure5.png','results/plots/supplementary_figure5a.png',
                #'results/plots/supplementary_figure5b.png','results/plots/supplementary_table1.tsv',
                #'results/plots/supplementary_table2.tsv','results/plots/supplementary_table3.tsv',
                #'results/plots/supplementary_figure4a.png','results/plots/supplementary_figure4b.png'
                 #lambda wildcards: expand("results/plots/hsa_PRJNA512378_{miRNA}_HeLa_noise_analysis.tsv", miRNA=config['PRJNA512378']['HeLa']['all_mirnas']),
                 #lambda wildcards: expand("results/plots/hsa_PRJNA231155_miR-137-3p_U251_noise_analysis.tsv"),
                 #lambda wildcards: expand("results/plots/hsa_PRJNA231155_miR-137-3p_U343_noise_analysis.tsv"),
                 #lambda wildcards: expand("results/plots/hsa_PRJNA292016_miR-141-3p_Du145_noise_analysis.tsv"),
                 #lambda wildcards: expand("results/plots/hsa_PRJNA304643_miR-1343-3p_A549_noise_analysis.tsv"),
                 #lambda wildcards: expand("results/plots/hsa_PRJNA304643_miR-1343-3p_16HBE14o_noise_analysis.tsv"),
                 #lambda wildcards: expand("results/plots/hsa_PRJNA223608_{miRNA}_U20S_noise_analysis.tsv", miRNA=config['PRJNA223608']['miRNAs']),
                 #lambda wildcards: expand("results/plots/mmu_PRJNA340017_{miRNA}_NMuMG_noise_analysis.tsv", miRNA=config['PRJNA340017']['miRNAs']),
                 #lambda wildcards: expand("results/plots/mmu_PRJNA309441_{miRNA}_CD4_noise_analysis.tsv", miRNA=config['PRJNA309441']['miRNAs']),
                 #lambda wildcards: expand("results/plots/mmu_PRJNA270999_miR-294-3p_ESCs_noise_analysis.tsv"),
                 'results/kallisto/SRR4054984/abundance.tsv','results/kallisto/SRR4054985/abundance.tsv',
                 'results/kallisto/SRR4055002/abundance.tsv','results/kallisto/SRR4055005/abundance.tsv'

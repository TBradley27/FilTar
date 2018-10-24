import sys
import os

if len(snakemake.input) == 1:
        os.system('salmon quant -p {} -i results/hsa_GRCh38_index/ -l A -r {} -o results/salmon/{}'.format(snakemake.threads,snakemake.input[0],snakemake.wildcards['accession']))
else:
        os.system('salmon quant -p {} -i results/hsa_GRCh38_index/ -l A -1 {} -2 {} -o results/salmon/{}'.format(snakemake.threads,snakemake.input[0],snakemake.input[1],snakemake.wildcards['accession']))                                                                                                                                                  

import sys
import os

if len(snakemake.input) == 2:
        os.system('salmon quant -p {} -i {} -l A -r {} -o results/salmon/{}'.format(snakemake.threads,snakemake.input[1],snakemake.input[0],snakemake.wildcards['accession']))
else:
        os.system('salmon quant -p {} -i {} -l A -1 {} -2 {} -o results/salmon/{}'.format(snakemake.threads,snakemake.input[2],snakemake.input[0],snakemake.input[1],snakemake.wildcards['accession']))                                                                                                                                                  

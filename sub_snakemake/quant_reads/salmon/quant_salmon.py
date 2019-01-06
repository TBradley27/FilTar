import sys
import os

index = snakemake.input['index'].replace('sa.bin','') # remove file subtring to leave just the path to the index 

print(index)

if len(snakemake.input) == 2:
        os.system('salmon quant -p {} -i {} -l A -r {} -o results/salmon/runs/{}'.format(snakemake.threads,index,snakemake.input['reads'][0],snakemake.wildcards['accession']))
else:
        os.system('salmon quant -p {} -i {} -l A -1 {} -2 {} -o results/salmon/runs/{}'.format(snakemake.threads,index,snakemake.input['reads'][0],snakemake.input['reads'][1],snakemake.wildcards['accession']))                                                                                                                                                  

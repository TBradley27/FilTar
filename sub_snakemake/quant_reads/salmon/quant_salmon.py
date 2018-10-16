import sys
import os

if len(snakemake.input) == 1:
        os.system('salmon quant -p {} -i {} -l A -r {} -o results/{}/'.format(snakemake.params['threads'],snakemake.input['index'],snakemake.index['reads'],snakemake.wildcards['accession']))
else:
        os.system('salmon quant -p {} -i {} -l A -1 {} -2 {} -o results/{}/'.format(snakemake.params['threads'],snakemake.input['index'],snakemake.index['forward_reads'],snakemake.index['reverse_reads'],snakemake.wildcards['accession']))
~                                                                                                                                                  

import sys
import os

if len(snakemake.input) == 2:
        os.system('salmon quant --validateMappings --rangeFactorizationBins 4 --seqBias --posBias --gcBuas -p {} -i {} -l A -r {} -o {}'.format(snakemake.threads,snakemake.input['index'],snakemake.input['reads'][0],snakemake.wildcards['species'],snakemake.wildcards['run_accession'],snakemake.output))
else:
        os.system('salmon quant --validateMappings --rangeFactorizationBins 4 --seqBias --posBias --gcBias -p {} -i {} -l A -1 {} -2 {} -o {}'.format(snakemake.threads,snakemake.input['index'],snakemake.input['reads'][0],snakemake.input['reads'][1],snakemake.wildcards['species'],snakemake.wildcards['run_accession'], snakemake.output))                                                                                                                                                  

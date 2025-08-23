import sys
import os

index_base_name = snakemake.input['index'].replace('.1.ht2','')

if len(snakemake.input['reads']) == 1:
	print('single_end_processing')
	os.system('hisat2 -x {} -p {} -U {} -S {} {}'.format(index_base_name,snakemake.threads,snakemake.input['reads'],snakemake.output,snakemake.params[0]))
else:
	print('paired_end_processing')
	os.system('hisat2 -x {} -p {} -1 {} -2 {} -S {} {}'.format(index_base_name,snakemake.threads,snakemake.input['reads'][0], snakemake.input['reads'][1], snakemake.output, snakemake.params[0]))

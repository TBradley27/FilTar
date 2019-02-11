import sys
import os

print(snakemake.threads)
print(snakemake.params[0])
print(snakemake.input[0])

index_base_name = snakemake.input['index'].replace('.1.ht2','')

if len(snakemake.input['reads']) == 1:
	print('single_end')
	os.system('hisat2 -x {} -p {} -U {} -S {} {} --summary-file reports/hisat2/{}.txt --new-summary'.format(index_base_name,snakemake.threads,snakemake.input['reads'],snakemake.output['sam_file'],snakemake.params[0], snakemake.wildcards['accession'],snakemake.wildcards['accession']))
else:
	print('paired_end')
	os.system('hisat2 -x {} -p {} -1 {} -2 {} -S {} {} --summary-file reports/hisat2/{}.txt --new-summary'.format(index_base_name,snakemake.threads,snakemake.input['reads'][0], snakemake.input['reads'][1], snakemake.output['sam_file'], snakemake.params[0], snakemake.wildcards['accession'],snakemake.wildcards['accession']))

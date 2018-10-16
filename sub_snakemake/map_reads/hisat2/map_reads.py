import sys
import os

print(snakemake.params[0])
print(snakemake.input[0])

if len(snakemake.input) == 1:
	os.system('hisat2 -x data/hsa -U {} -S {} {}'.format(snakemake.input,snakemake.output,snakemake.params[0]))
else:
	os.system('hisat2 -x data/hsa -1 {} -2 {} -S {} {}'.format(snakemake.input[0], snakemake.input[1], snakemake.output, snakemake.params[0]))

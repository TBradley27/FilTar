import os

if len(snakemake.input) > 1:
	os.system("bedtools unionbedg -i {} > {}".format(snakemake.input, snakemake.output))
else:
	os.system("cp {} {}".format(snakemake.input, snakemake.output))

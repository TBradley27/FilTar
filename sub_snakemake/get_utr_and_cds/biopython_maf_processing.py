#!/bin/env python

import sys
import re

# replace "./alignio-maf" with the full path of the alignio-maf branch
# you cloned from github, or keep if it's in the current directory

#sys.path.insert(1, "/Users/bradleyt/Documents/testing/alignio-maf")

from Bio import AlignIO
from Bio.AlignIO import MafIO

print(snakemake.output[0])
print(snakemake.input[0])
print(snakemake.wildcards['species'])
print(snakemake.wildcards['chrom'])


if snakemake.wildcards['species'] == "hsa":
       build = "hg38" 
elif snakemake.wildcards['species'] == "mmu":
	build = "mm10"
else:
	build = ''

print(build)

idx = AlignIO.MafIO.MafIndex(snakemake.output[0],  snakemake.input[0], "{}.chr{}".format(build, snakemake.wildcards['chrom'])  )

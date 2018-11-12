#!/bin/env python

import sys
import re

# replace "./alignio-maf" with the full path of the alignio-maf branch
# you cloned from github, or keep if it's in the current directory

#sys.path.insert(1, "/Users/bradleyt/Documents/testing/alignio-maf")

from Bio import AlignIO
from Bio.AlignIO import MafIO

if snakemake.wildcards['species'] == "hsa":
       build = "hg38" 
elif snakemake.wildcards['species'] == "mmu":
	build = "mm10"
else:
	build = ''

print(sys.argv)

idx = AlignIO.MafIO.MafIndex(snakemake.output[0],  snakemake.input[1], "{}.chr{}".format(build, snakemake.wildcards['chrom'])  )

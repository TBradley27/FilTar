#!/bin/env python3

import sys
import re
from Bio import AlignIO
from Bio.AlignIO import MafIO

if snakemake.wildcards['species'] == "hsa":
       build = "hg38" 
elif snakemake.wildcards['species'] == "mmu":
	build = "mm10"
else:
	build = ''

idx = AlignIO.MafIO.MafIndex(snakemake.output[0],  snakemake.input[0], "{}.chr{}".format(build, snakemake.wildcards['chrom'])  )

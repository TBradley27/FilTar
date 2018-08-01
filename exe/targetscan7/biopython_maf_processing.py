#!/bin/env python

import sys
import re

# replace "./alignio-maf" with the full path of the alignio-maf branch
# you cloned from github, or keep if it's in the current directory

#sys.path.insert(1, "/Users/bradleyt/Documents/testing/alignio-maf")

from Bio import AlignIO
from Bio.AlignIO import MafIO

if sys.argv[3] == "hsa":
       build = "hg38" 
elif sys.argv[3] == "mmu":
	build = "mm10"
else:
	build = ''

print(sys.argv)

idx = AlignIO.MafIO.MafIndex(sys.argv[1],  sys.argv[2], "{}.chr{}".format(build, sys.argv[4])  )

#!/bin/env python

#    FilTar: Integrating RNA-Seq data to improve microRNA target prediction accuracy in animals
#    Copyright (C) 2019 Thomas Bradley
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>

import sys
import functools
from functools import reduce
import re
import os
from Bio import AlignIO
from Bio.AlignIO import MafIO
from subprocess import call

def add_alignment():
           global start_pos
           global end_pos
           if strand == -1:
               start_pos = start_pos[::-1] # Reverse Order
               end_pos = end_pos[::-1]
           else:
               pass

           new_multiple_alignment = idx.get_spliced(start_pos, end_pos, strand) # splice through the index
           for i in range(len(new_multiple_alignment)):
              new_multiple_alignment[i].id = re.sub('\..*','', new_multiple_alignment[i].id) # strip chomosome information from id
              new_multiple_alignment[i].id = accession + '\t' + new_multiple_alignment[i].id # label all alignments with relevant transcript accession
           
           global big_alignment   
           big_alignment += new_multiple_alignment.format('fasta')
           return()

if os.stat(snakemake.input['bed']).st_size == 0:
     target = open(snakemake.output[0], 'w')
else:
    if snakemake.wildcards['species'] == "hsa":    #Identify the species identifier passed through the command line
       build = "hg38"
    elif snakemake.wildcards['species'] == "mmu":
       build = "mm10"
    else:
       build = ''

    idx = AlignIO.MafIO.MafIndex(snakemake.input['maf_index'], snakemake.input['maf'], "{}.chr{}".format(build, snakemake.wildcards['chrom'])  )

    start_pos = []
    end_pos = []
    accession = 'empty'
    with open(snakemake.input['bed'] ) as f:
       big_alignment = ''
       for line in f:    #Open and loop through line-by-line the relevant BED file
          parts = line.split()  

          if accession == 'empty':
             accession = parts[4]   
             start_pos.append(int(parts[1]))
             end_pos.append(int(parts[2]))
             strand = (int(parts[3]))

          elif parts[4] == accession:       # If accession has multiple entries in bed, add additional genome co-ordinates to rel. lists
             start_pos.append(int(parts[1]))
             end_pos.append(int(parts[2]))
             strand = (int(parts[3]))

          else:
             add_alignment()
   
             start_pos = [] # initialise a new transcript record
             end_pos = []
             start_pos.append(int(parts[1]))
             end_pos.append(int(parts[2]))
             strand = (int(parts[3]))
             accession = parts[4]
       else: #Not sure, but I think this is for transcripts with only one 'exon', which are not the first transcript in the bed record.

         add_alignment()     

       result = reduce(lambda x, y: x.replace(y, snakemake.config["TaxID"][y]), snakemake.config["TaxID"], big_alignment) # get NCBI taxonomic IDs
       result = result.replace('\n','') # convert from fasta to tsv
       result = result.replace('>','\n') # convert from fasta to tsv
       result = re.sub('(\s[0-9]{4,7})',r'\1\t',result) # convert from fasta to tsv
       result = re.sub('\n','',result, count=1) # remove leading empty line

       target = open(snakemake.output[0], 'w')

       for line in iter(result.splitlines()):
          pattern = re.compile('\.[0-9][0-9]?\sN+')
          pattern2 = re.compile('T[0-9]+\sN+')
          pattern3 = re.compile('^N+$') # when ref transcript is unknown - it leaves a trailing lines of Ns which need to be removed
          if 'delete' in line:
              pass
          elif 'unknown' in line:
              pass
          elif pattern.search(line):
              pass
          elif pattern2.search(line):
              pass
          elif pattern3.search(line):
              pass
          else:
              line2 = line + '\n'
              target.write(line2)

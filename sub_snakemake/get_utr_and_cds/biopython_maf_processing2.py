#!/bin/env python

import sys
import functools
from functools import reduce
import re
import os
from Bio import AlignIO
from Bio.AlignIO import MafIO
from subprocess import call

print(snakemake.config["TaxID"]['hg38'])

if snakemake.wildcards['species'] == "hsa":    #Identify the species identifier passed through the command line
       build = "hg38"
elif snakemake.wildcards['species'] == "mmu":
        build = "mm10"
else:
        build = ''

print (snakemake)
print (build)

idx = AlignIO.MafIO.MafIndex(snakemake.input['maf_index'], snakemake.input['maf'], "{}.chr{}".format(build, snakemake.wildcards['chrom'])  )

start_pos = []
end_pos = []
accession = 'empty'
with open(snakemake.input['bed'] ) as f:
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
           # write process - triggered once script hits a record != accession

           if strand == -1:
               start_pos = start_pos[::-1] # Reverse Order
               end_pos = end_pos[::-1]
           else:
               pass

           print (accession)
           print (start_pos)
           print (end_pos)
           print (strand)             

           new_multiple_alignment = idx.get_spliced(start_pos, end_pos, strand) # splice through the index
           AlignIO.write(new_multiple_alignment, "results/{}.fa".format(accession), "fasta")
           
           call("exe/targetscan7/convert_fasta_to_tsv.sh results/{}.fa {} > tmp{}.tsv".format(accession, accession, snakemake.wildcards['chrom']), shell=True) #Execute a shell command
           call(["rm", "results/{}.fa".format(accession)])
           call("cat tmp{}.tsv >> results/hsa_chr{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['chrom']), shell=True)
           call(["rm","tmp{}.tsv".format(snakemake.wildcards['chrom'])])
           
           start_pos = [] # initialise a new transcript record
           end_pos = []
           start_pos.append(int(parts[1]))
           end_pos.append(int(parts[2]))
           strand = (int(parts[3]))
           accession = parts[4]
    else: #Not sure, but I think this is for transcripts with only one 'exon', which are not the first transcript in the bed record.

       if strand == -1:
               start_pos = start_pos[::-1] # Reverse Order
               end_pos = end_pos[::-1]
       else:        
               pass

       print (accession)
       print (start_pos)
       print (end_pos)
       new_multiple_alignment = idx.get_spliced(start_pos, end_pos, strand)
       AlignIO.write(new_multiple_alignment, "results/{}.fa".format(accession), "fasta")

       call("exe/targetscan7/convert_fasta_to_tsv.sh results/{}.fa {} > tmp{}.tsv".format(accession, accession, snakemake.wildcards['chrom']), shell=True) #Execute a shell command
       call(["rm", "results/{}.fa".format(accession)])
       call("cat tmp{}.tsv >> results/hsa_chr{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['chrom']), shell=True)
       call(["rm","tmp{}.tsv".format(snakemake.wildcards['chrom'])])

f = open("results/hsa_chr{}_msa_tmp.tsv".format(snakemake.wildcards['chrom']), 'r')

# prevents writing to an already existing file
call(["rm",snakemake.output[0]])

target = open(snakemake.output[0], 'w')

for line in iter(f):
	result = reduce(lambda x, y: x.replace(y, snakemake.config["TaxID"][y]), snakemake.config["TaxID"], line)
	pattern = re.compile('\.[0-9][0-9]?\sN+')
	if 'delete' in result: # deletes some species to get an 84-way alingment instead of a 100-way alignment
		pass
	elif 'unknown' in result:
		pass
	elif pattern.search(result): # deletes malprocessed lines in which repeat Ns comprise the second column
		pass
	else:
		target.write(result)
	del result
f.close()

call(["rm","results/hsa_chr{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'])])

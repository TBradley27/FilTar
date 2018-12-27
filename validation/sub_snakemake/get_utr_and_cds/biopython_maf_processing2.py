#!/bin/env python

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
           # write process - triggered once script hits a record != accession

          #AlignIO.write(new_multiple_alignment, "results/{}_{}.fa".format(accession, snakemake.wildcards['feature']), "fasta")
           
<<<<<<< HEAD:validation/sub_snakemake/get_utr_and_cds/biopython_maf_processing2.py
           call("exe/targetscan7/convert_fasta_to_tsv.sh results/{}_{}.fa {} > tmp{}_{}.tsv".format(accession,snakemake.wildcards['feature'],  accession,snakemake.wildcards['chrom'], snakemake.wildcards['feature']), shell=True) #Execute a shell command
           call(["rm", "results/{}_{}.fa".format(accession, snakemake.wildcards['feature'])])
           call("cat tmp{}_{}.tsv >> results/hsa_chr{}_{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'], snakemake.wildcards['chrom'], snakemake.wildcards['feature']), shell=True)
           call(["rm","tmp{}_{}.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'])])
=======
           #call("exe/targetscan7/convert_fasta_to_tsv.sh results/{}_{}.fa {} > tmp{}_{}.tsv".format(accession, snakemake.wildcards['feature'], accession, snakemake.wildcards['chrom'], snakemake.wildcards['feature']), shell=True) #Execute a shell command
           #call(["rm", "results/{}_{}.fa".format(accession, snakemake.wildcards['feature'])])
           #call("cat tmp{}_{}.tsv >> results/hsa_chr{}_{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'], snakemake.wildcards['chrom'], snakemake.wildcards['feature']), shell=True)
           #call(["rm","tmp{}_{}.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'])])
>>>>>>> c6b51b5... refactor code for generating alignments:sub_snakemake/get_utr_and_cds/biopython_maf_processing2.py
           
           start_pos = [] # initialise a new transcript record
           end_pos = []
           start_pos.append(int(parts[1]))
           end_pos.append(int(parts[2]))
           strand = (int(parts[3]))
           accession = parts[4]
    else: #Not sure, but I think this is for transcripts with only one 'exon', which are not the first transcript in the bed record.

      add_alignment()     

       #call("exe/targetscan7/convert_fasta_to_tsv.sh results/{}_{}.fa {} > tmp{}_{}.tsv".format(accession, snakemake.wildcards['feature'], accession, snakemake.wildcards['chrom'], snakemake.wildcards['feature']), shell=True) #Execute a shell command
       #call(["rm", "results/{}_{}.fa".format(accession, snakemake.wildcards['feature'])])
       #call("cat tmp{}_{}.tsv >> results/hsa_chr{}_{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'], snakemake.wildcards['chrom'], snakemake.wildcards['feature']), shell=True)
       #call(["rm","tmp{}_{}.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'])])

#f = open("results/hsa_chr{}_{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature']), 'r')

# prevents writing to an already existing file
#call(["rm",snakemake.output[0]])

result = reduce(lambda x, y: x.replace(y, snakemake.config["TaxID"][y]), snakemake.config["TaxID"], big_alignment) # get NCBI taxonomic IDs
result = result.replace('\n','') # convert from fasta to tsv
result = result.replace('>','\n') # convert from fasta to tsv
result = re.sub('(\s[0-9]{4,6})',r'\1\t',result) # convert from fasta to tsv
result = re.sub('\n','',result, count=1) # remove leading empty line

target = open(snakemake.output[0], 'w')

for line in iter(result.splitlines()):
        pattern = re.compile('\.[0-9][0-9]?\sN+')
        pattern2 = re.compile('T[0-9]+\sN+')
        if 'delete' in line:
            pass
        elif 'unknown' in line:
            pass
        elif pattern.search(line):
            pass
        elif pattern2.search(line):
            pass
        else:
            line2 = line + '\n'
            target.write(line2)

##for line in iter(f):
#	result = reduce(lambda x, y: x.replace(y, snakemake.config["TaxID"][y]), snakemake.config["TaxID"], line)
#	pattern = re.compile('\.[0-9][0-9]?\sN+')
#	pattern2 = re.compile('T[0-9]+\sN+') # in cases in which the transcript identifer does not have a version number
#	if 'delete' in result: # deletes some species to get an 84-way alingment instead of a 100-way alignment
#		pass
#	elif 'unknown' in result:
#		pass
#	elif pattern.search(result): # deletes malprocessed lines in which repeat Ns comprise the second column
#		pass
#	elif pattern2.search(result): 
#		pass
#	else:
#		target.write(result)
#	del result
#f.close()

#call(["rm","results/hsa_chr{}_{}_msa_tmp.tsv".format(snakemake.wildcards['chrom'], snakemake.wildcards['feature'])])

import sys
import os

if len(snakemake.input['fastq']) == 1:
        print('single-end processing')
        os.system('kallisto quant --bias -i {} -t {} -o results/kallisto/{}/ --single -l {} -s {} {} 2> {}'.format(snakemake.input['index'],snakemake.wildcards['accession'],snakemake.params['frag_length_mean'],snakemake.params['frag_length_sd'],snakemake.input['fastq'][0],snakemake.log[0]))
else:
        print('paired-end processing')
        os.system('kallisto quant --bias -i {} -t {} -o results/kallisto/{}/ {} {} 2> {}'.format(snakemake.input['index'],snakemake.wildcards['accession'],snakemake.input['fastq'][0],snakemake.input['fastq'][1],snakemake.log[0]))                                                                                          

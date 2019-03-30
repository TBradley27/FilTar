import sys
import os

if 'single_end' in snakemake.output[0]:
	os.system("wget -nv --directory-prefix=data/single_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/00{}/{}/{}.fastq.gz || wget -nv --directory-prefix=data/single_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{}/{}.fastq.gz".format(snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'][-1],snakemake.wildcards['accession'],snakemake.wildcards['accession'],snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'],snakemake.wildcards['accession']))
elif 'paired_end' in snakemake.output[0]:
	os.system("wget -nv --directory-prefix=data/paired_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/00{}/{}/{}_{}.fastq.gz || wget -nv --directory-prefix=data/paired_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{}/{}_{}.fastq.gz".format(snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'][-1],snakemake.wildcards['accession'],snakemake.wildcards['accession'],snakemake.wildcards['mate_number'],snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'],snakemake.wildcards['accession'],snakemake.wildcards['mate_number']))

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
import os

if 'single_end' in snakemake.output[0]:
	os.system("wget -nv --directory-prefix=data/single_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/00{}/{}/{}.fastq.gz || wget -nv --directory-prefix=data/single_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{}/{}.fastq.gz".format(snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'][-1],snakemake.wildcards['accession'],snakemake.wildcards['accession'],snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'],snakemake.wildcards['accession']))
elif 'paired_end' in snakemake.output[0]:
	os.system("wget -nv --directory-prefix=data/paired_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/00{}/{}/{}_{}.fastq.gz || wget -nv --directory-prefix=data/paired_end/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{}/{}_{}.fastq.gz".format(snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'][-1],snakemake.wildcards['accession'],snakemake.wildcards['accession'],snakemake.wildcards['mate_number'],snakemake.wildcards['accession'][0:6],snakemake.wildcards['accession'],snakemake.wildcards['accession'],snakemake.wildcards['mate_number']))

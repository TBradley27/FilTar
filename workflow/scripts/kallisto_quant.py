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

if len(snakemake.input['fastq']) == 1:
        print('single-end processing')
        os.system('kallisto quant --bias -i {} -t {} -o results/kallisto/{}/ --single -l {} -s {} {} 2> {}'.format(snakemake.input['index'],snakemake.wildcards['accession'],snakemake.params['frag_length_mean'],snakemake.params['frag_length_sd'],snakemake.input['fastq'][0],snakemake.log[0]))
else:
        print('paired-end processing')
        os.system('kallisto quant --bias -i {} -t {} -o results/kallisto/{}/ {} {} 2> {}'.format(snakemake.input['index'],snakemake.wildcards['accession'],snakemake.input['fastq'][0],snakemake.input['fastq'][1],snakemake.log[0]))                                                                                          

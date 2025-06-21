#FilTar: Integrating RNA-Seq data to improve microRNA target prediction accuracy in animals
#Copyright (C) 2019 Thomas Bradley
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

import re
import pandas

# import metadata
metadata = pandas.read_table("metadata.tsv")

configfile: "config/basic.yaml"
configfile: "config/species.yaml"
configfile: "config/dependencies.yaml"

# Include modules with rules that other modules depend on first
include: "modules/data_download/Snakefile"
include: "modules/trim_reads/Snakefile"
include: "modules/quant_reads/salmon/Snakefile"
include: "modules/mirna/Snakefile"
include: "modules/get_target_coordinates/Snakefile"

if config['reannotation'] == True:
    include: "modules/with_reannotation/Snakefile"
    include: "modules/with_reannotation/map_reads/hisat2/Snakefile"
elif config['reannotation'] == False:
    include: "modules/without_reannotation/Snakefile"
else:
    raise Exception("\nPlease enter a value of either 'True' or 'False' for the 'reannotation' key. Default values can be set in config/basic.yaml\n")

if config['conservation'] == True:
    include: "modules/get_utr_and_cds/with_conservation/Snakefile"
elif config['conservation'] == False:
    include: "modules/get_utr_and_cds/without_conservation/Snakefile"
else:
    raise Exception("\nPlease enter a value of either 'True' or 'False' for the 'conservation' key. Default values can be set in config/basic.yaml\n")

if config['sequence_data_source'] == 'ENA':
    include: "modules/data_download/ENA/Snakefile"
elif config['sequence_data_source'] == 'SRA':
    include: "modules/data_download/SRAtoolkit/Snakefile"
elif config['sequence_data_source'] == 'User':
    pass
else:
    raise Exception("\nPlease enter a value of either 'ENA' or 'SRA' or 'User' for the 'sequence_data_source' key. Default values can be set in config/basic.yaml\n")

# Include target prediction algorithm Snakefiles last
if config['prediction_algorithm'] == 'TargetScan7':
    pass
elif config['prediction_algorithm'] == 'miRanda':
    pass
else:
    raise Exception("\nPlease enter a valid name for a miRNA target prediction algorithm. Choose either 'TargetScan7' or 'miRanda'\n")

if config['prediction_algorithm'] == 'TargetScan7' and config['reannotation'] == True:
    include: "modules/target_prediction/targetscan/Snakefile"
    include: "modules/target_prediction/targetscan/with_reannotation/Snakefile"
elif config['prediction_algorithm'] == 'TargetScan7' and config['reannotation'] == False:
    include: "modules/target_prediction/targetscan/Snakefile"
    include: "modules/target_prediction/targetscan/without_reannotation/Snakefile"
elif config['prediction_algorithm'] == 'miRanda' and config['reannotation'] == True:
    include: "modules/target_prediction/miRanda/Snakefile"
    include: "modules/target_prediction/miRanda/with_reannotation/Snakefile"
elif config['prediction_algorithm'] == 'miRanda' and config['reannotation'] == False:
    include: "modules/target_prediction/miRanda/Snakefile"
    include: "modules/target_prediction/miRanda/without_reannotation/Snakefile"

if config['conservation'] == True and config['prediction_algorithm'] == 'miRanda':
    raise Exception("miRanda cannot be used when the conservation option is set to True")
else:
    pass

if config['reannotation'] == True and config['context'] == 'reference':
    raise Exception("The reannotation option cannot be set to True and context option set to 'reference' at the same time")
else:
    pass

for transcript in list(config['transcripts']):
    if re.match('^ENS[A-Z]+[0-9]+.[1-9]{1,2}$',transcript):
        pass
    else:
        raise Exception('\nInvalid transcript identifier "{}". Identifiers must adhere to official Ensembl identifier patterns e.g. "ENSMUST00000189888.6". Please revise.\n'.format(transcript))

wildcard_constraints:
    species="[a-z]{3,4}",
    tissue_tx_model="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue_tx_model wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    tissue_tx_exp="((?!chr([A-Z]|\d)).)*", # pattern to ensure tissue_tx_exp wildcard does not contain the following pattern: chr[0-9] or chr[A-Z]
    chrom="[A-Za-z0-9]{1,5}",
    feature="(3UTR|CDS)",
    ensembl_release="[0-9]{2,3}",
    genus_species="[A-Z][a-z]+_[a-z]+"

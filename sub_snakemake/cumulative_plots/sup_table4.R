library(tidyverse)

### read in the data

U251 = read_tsv(snakemake@input['U251'])
U343 = read_tsv(snakemake@input['U343'])
DU145 = read_tsv(snakemake@input['Du145'])
A549 = read_tsv(snakemake@input['A549'])
16HBE14o = read_tsv(snakemake@input['16HBE14o'])
HeLa = read_tsv(snakemake@input['HeLa'])
U20S = read_tsv(snakemake@input['U20S'])
Kidney = read_tsv(snakemake@input['Kidney'])
Lung = read_tsv(snakemake@input['Lung'])
Skeletal_muscle = read_tsv(snakemake@input['Skeletal_muscle'])
Thyroid = read_tsv(snakemake@input['Thyroid'])
Bone_marrow = read_tsv(snakemake@input['Bone_marrow'])
NMuMG = read_tsv(snakemake@input['NMuMG'])
ESCs = read_tsv(snakemake@input['ESCs'])
CD4 = read_tsv(snakemake@input['CD4'])

library(plyr)
library(tidyverse)
source('validation/sub_snakemake/canonical_targets/get_canonical_AIR_function.R')

AIRs = get_canonical_AIRs(snakemake@input[[1]])

write.table(AIRs, snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

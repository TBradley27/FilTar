library(plyr)
library(tidyverse)
source('fix_ts_output_function.R')

ts_sites = fix_ts_sites(snakemake@input[[1]])


write.table(ts_sites,
file=snakemake@output[[1]],
sep="\t",
row.names=FALSE,
quote=FALSE)

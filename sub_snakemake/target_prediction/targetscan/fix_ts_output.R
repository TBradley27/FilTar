library(plyr)
library(tidyverse)
library(filtar)

ts_sites = fix_ts_output(snakemake@input[[1]])

write.table(ts_sites,
file=snakemake@output[[1]],
sep="\t",
row.names=FALSE,
quote=FALSE)

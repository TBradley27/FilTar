library(readr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

three_utrs = read_tsv(args[1], col_names=c('chrom','start','stop','transcript_id','dummy','sign'))
CDS = read_tsv(args[2], col_names=c('chrom','start','stop','transcript_id','dummy', 'sign'))

CDS = CDS[CDS$transcript_id %in% three_utrs$transcript_id,]

write.table(CDS, args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

#!/bin/env Rscript

united_quant = readr::read_tsv(
        file=snakemake@input[[1]],
        col_types=cols(.default = 'd', Name = 'c')
)

united_quant2 = filtar::AvgSalmonQuant(united_quant)

avg_quant = united_quant2[,c('Name','avg')]
colnames(avg_quant) = c('Name','TPM')

write.table(
        x=avg_quant,
        file=snakemake@output[[1]],
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE,
        sep="\t"
)      

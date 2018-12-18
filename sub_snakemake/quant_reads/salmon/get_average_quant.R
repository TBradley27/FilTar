#!/bin/env Rscript

united_quant = readr::read_tsv(
        file=snakemake@input[[1]],
        col_types=readr::cols(.default = 'd', Name = 'c')
)

united_quant2 = filtar::AvgSalmonQuant(united_quant)

avg_quant = tibble::tibble(Name=united_quant2$Name, Length=20, EffectiveLength=20.00, TPM=united_quant2$avg, NumReads=20.00)

real_output = paste(snakemake@output[[1]],'quant.sf',sep="")

write.table(
        x=avg_quant,
        file=real_output,
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE,
        sep="\t"
)      

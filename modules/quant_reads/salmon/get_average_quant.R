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

united_quant = readr::read_tsv(
        file=snakemake@input[[1]],
        col_types=readr::cols(.default = 'd', Name = 'c')
)

united_quant2 = filtar::AvgSalmonQuant(united_quant)

avg_quant = tibble::tibble(Name=united_quant2$Name, Length=20, EffectiveLength=20.00, TPM=united_quant2$avg, NumReads=20.00)

real_output = paste(snakemake@output[[1]],'quant.sf',sep="/")
dir.create(snakemake@output[[1]])

write.table(
        x=avg_quant,
        file=real_output,
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE,
        sep="\t"
)      

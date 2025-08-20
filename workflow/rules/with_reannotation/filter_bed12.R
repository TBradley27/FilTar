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

input = readr::read_tsv(snakemake@input[[1]], col_names=FALSE, col_types='ciiciciiiicc')

if (length(snakemake@config[['transcripts']]) == 0) {
        file.copy(from=snakemake@input[[1]],to=snakemake@output[[1]])
} else {
	transcripts = gsub('\\..*','', snakemake@config[['transcripts']])

	filtered_input = input[input$X4 %in% transcripts,]

	write.table(filtered_input,snakemake@output[[1]], col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}

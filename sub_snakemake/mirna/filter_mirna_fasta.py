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

from Bio import SeqIO

records = list(SeqIO.parse(snakemake.input[0],'fasta'))
filtered_records = []

for record in records:
	if record.id in snakemake.config['mirnas']:
		filtered_records.append(record)
	elif len(snakemake.config['mirnas']) == 0: # use all records if mirna config entry is empty
		filtered_records.append(record)
	else:
		pass

SeqIO.write(filtered_records,snakemake.output[0],'fasta')

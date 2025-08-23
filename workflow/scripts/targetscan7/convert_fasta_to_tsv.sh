#!/bin/bash

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

# Convert text files from fasta format to TargetScan7 format
err() {
  echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: $@" >&2; exit 1
}

fasta_file="$1"
accession=$2

cat "$fasta_file"                             |
sed 's/\./ /1'                                | # Separate species from chr num
sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g'    | # Remove newline characters
sed 's/>//1'                                  | # Not sure of function, if any
sed 's/>/\n/g'                                | # Ensure one line per sequence
sed 's/[[:space:]]//3g'                       | # Remove spaces in sequence
awk -v accession="$2" '{ OFS="\t" }{ print accession,$1,$3 }' 

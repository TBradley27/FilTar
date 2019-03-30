#!/bin/bash

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

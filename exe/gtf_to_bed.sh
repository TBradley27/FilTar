#!/bin/bash

# Convert a file from GTF format to BED format for use by either BedTools
# or biopython. 

# TODO(bradleyt): Write repeated sed commands as a function
# TODO(bradleyt): Ensure path variables are not too long for single line

gtf="$1"; output="$2" # aliases

grep three_prime_utr "$gtf" |
awk '{print $1,$4,$5,$7,$14$16}' |
  sed 's/+/1/g'   |
  sed 's/-/-1/g'  |
  sed 's/\"//g'   |
  sed 's/;/./g'   |
  sed 's/\.//2'   |
  tr ' '  \\t     |
  awk '{ OFS="\t" }{print $1,$2-1,$3,$4,$5}' > "$output"
    
  if [[ "$?" -ne 0 ]]; then
    err "Unable to complete biopython formatting" >&2
  fi

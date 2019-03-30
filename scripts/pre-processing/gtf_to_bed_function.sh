#!/bin/bash

gtf_to_bed () {

	input="$1"; tx_feature="$2"; output="$3"

	if [[ "$tx_feature" == '3UTR' ]]; then
		tx_feature_long=three_prime_utr
	elif [[ "$tx_feature" == 'CDS' ]]; then
		tx_feature_long=CDS
	else
	  echo "$2 is not a valid genomic feature"; exit 1
	fi

	grep "${tx_feature_long}" "$input" |
	grep -E "\s${tx_feature_long}\s"|         # Braces needed for correct search - test for white space either side of the pattern
	sed 's/^chr//g' |
	awk '{print $1,$4,$5,$7,$14$16}' |
	sed 's/+/1/g'   |
	sed 's/-/-1/g'  |
	sed 's/\"//g'   |
	sed 's/;/./g'   |
	sed 's/\.//2'   |
	tr ' '  \\t     |
	awk '{ OFS="\t" }{print $1,$2-1,$3,$4,$5}'   > "$output"
	    
	if [[ "$?" -ne 0 ]]; then
	   err "Unable to complete biopython formatting" >&2
	fi
}



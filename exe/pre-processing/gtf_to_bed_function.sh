#!/bin/bash

gtf_to_bed () {

	input="$1"; tx_feature="$2"; chromosome="$3"; output="$4"

	if [[ "$tx_feature" == '3UTR' ]]; then
		tx_feature_long=three_prime_utr
	elif [[ "$tx_feature" == 'CDS' ]]; then
		tx_feature_long=CDS
	else
	  echo "$2 is not a valid genomic feature"; exit 1
	fi

	grep "${tx_feature_long}" "$input" |
	grep -E "\s${tx_feature_long}\s" |         # Braces needed for correct search - test for white space either side of the pattern
	sed 's/^chr//g' |
	grep '^'$chromosome'\s' > tmp"$chromosome"  # CCDS match confounds CDS search - warning: mmu and hsa have different prefixes for this

	if [[ "$?" -ne 0 ]] || [[ ! -s tmp"$chromosome"  ]]; then
	  echo "Unable to complete initial processing" >&2; exit 1
	fi

	awk '{print $1,$4,$5,$7,$14$16}' tmp"$chromosome" |
	sed 's/+/1/g'   |
	sed 's/-/-1/g'  |
	sed 's/\"//g'   |
	sed 's/;/./g'   |
	sed 's/\.//2'   |
	tr ' '  \\t     |
	#awk '$2!=$3'    |  # the next line converts from 1 to 0-based indexing
	awk '{ OFS="\t" }{print $1,$2-1,$3,$4,$5}'   > "$output" # $6
	    
	if [[ "$?" -ne 0 ]]; then
	   err "Unable to complete biopython formatting" >&2
	fi
	  
	rm tmp"$chromosome"

}



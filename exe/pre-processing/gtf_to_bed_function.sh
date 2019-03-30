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
#	grep '^'$chromosome'\s' |  # CCDS match confounds CDS search - warning: mmu and hsa have different prefixes for this
	awk '{print $1,$4,$5,$7,$14$16}' |
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
}



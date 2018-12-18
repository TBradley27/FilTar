#!/bin/bash

function prep_alt_miRNA () {

	if [ $2 == 'hsa' ]
		then tax_id='9606'
	elif [ $2 == 'mmu' ]
		then tax_id='10090'
	else
		tax_id='foo'
	fi

	awk '{print $1}' $1 |
	grep -A 1 "$2-"  |
	sed 's/--//g'    |
	sed '/^$/d'      |
	sed 'N;s/\n/ /'  |
	sed 's/>//g'     |
	awk -v tax_id="$tax_id" '{ OFS="\t" }{ print $1,tax_id,$1,$2}' |
	grep -E '\s[ACTGU]{18,24}$' | 
	sed 's/$2-//g'
}

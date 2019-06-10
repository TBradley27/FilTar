cat $1 | grep '>' | awk '{{ print $1,"Ensembl","94",$7 }}' | tr ' ' '\t' | sed 's/gene_symbol://g' | sed 's/>//g'

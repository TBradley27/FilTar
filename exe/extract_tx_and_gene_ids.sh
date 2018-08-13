cat $1 | grep '>' | awk '{ print $1,"ensembl",'92',$7 }' | sed 's/\..//g' | tr ' ' '\t' | sed 's/gene_symbol://g' | sed 's/>//g'

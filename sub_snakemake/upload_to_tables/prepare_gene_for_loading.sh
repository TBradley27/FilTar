cat $1 | awk '$3="gene"' | awk '{OFS="\t"}{print $20}' | sed 's/"//g' | sed 's/;//g' | grep -v '^[ensembl|havana]' | sed '/^$/d' | sort | uniq

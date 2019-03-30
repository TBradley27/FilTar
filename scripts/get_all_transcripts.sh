cat $1 | awk '$3 == "transcript"' | awk '{OFS="\t"}{print $1,$14,$16}' | sed 's/"//g' | sed 's/;//g' | sed 's/\t/\./2' | grep -E "^$2	" | awk '{print $2}'

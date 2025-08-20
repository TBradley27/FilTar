awk '{OFS="\t"}{ print $1,$2,$3,$5,$5,$4}' $1 | sed 's/\t1$/\t+/g' | sed 's/\t-1$/\t-/g' > $2

sed 's/(+)//g' $1 | sed 's/(-)//g' | awk -v species=$2 '{OFS="\t"}{ print $1,species,$2}' > $3

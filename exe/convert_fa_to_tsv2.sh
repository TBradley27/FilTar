sed 's/(+)//g' $1 | sed 's/(-)//g' | sed 's/(-)//g' | tr '\n' '\t' | sed 's/>/\n/g' | sed '/^$/d' | awk -v species=$2 '{OFS="\t"}{ print $1,species,$2}' > $3

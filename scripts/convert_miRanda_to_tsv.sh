grep 'hit' $1
command=$? # get exit code

if [[ $command -eq 0  ]]
then
	grep -A 1 'hit' $1 | sed '/--/d' | grep '>' | sed 's/>//' > $2
else
	touch $2
fi

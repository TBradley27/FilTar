cat $1 | awk '$3 == "transcript"' | awk '{print $14,$16}' | sed 's/"//g' | sed 's/;//g' | tr ' ' '.'

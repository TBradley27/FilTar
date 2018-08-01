#!/bin/bash

cat $1 | awk '$3="transcript"' | awk '{OFS="\t"}{print $14,"Ensembl",92,$10}' | grep '^"ENST' | sort | uniq | sed 's/"//g' | sed  's/;//g'

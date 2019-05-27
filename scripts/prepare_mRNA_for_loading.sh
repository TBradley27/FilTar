#!/bin/bash

cat $1 | awk '$3="transcript"' | awk '{OFS="\t"}{print $14,"Ensembl",94,$10}' | sed 's/"//g' | grep -P '^ENS[A-Z]{1,4}000' | sort | uniq | sed  's/;//g'

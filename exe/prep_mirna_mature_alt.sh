#!/bin/bash

awk '{print $1}' $1 |
grep -A 1 'hsa-' |
sed 's/--//g'    |
sed '/^$/d'      |
sed 'N;s/\n/ /'  |
sed 's/>//g'     |
awk '{ OFS="\t" }{ print $1,'9606',$1,$2}'  |
sed 's/hsa-//g'

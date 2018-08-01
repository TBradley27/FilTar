#!/bin/bash

sed -e '38 s;PCT_param;src\/targetscan7\/PCT_param;g' $1 | sed '41i if ($ARGV[1] eq 'hsa') {$refGenome = 9606;} else {$refGenome = 10090;}'

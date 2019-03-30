#!/bin/bash

sed -e '77 s;TA_SPS_by;exe\/targetscan7\/TA_SPS_by;g' $1 | sed -e '80 s;Agarwal_;exe\/targetscan7\/Agarwal_;g' | sed -e '83 s;All;exe\/targetscan7\/All;g' | sed -e '112 s;RNAplfold_in;exe\/targetscan7\/RNAplfold_in;g' |  sed '116i if ($ARGV[1] eq 'hsa') {$REF_SPECIES = 9606;} else {$REF_SPECIES = 10090;}'


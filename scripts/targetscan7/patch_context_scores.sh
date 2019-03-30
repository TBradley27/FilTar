#!/bin/bash

sed -e '80 s;TA_SPS_by;scripts\/targetscan7\/TA_SPS_by;g' $1 | sed -e '83 s;Agarwal_;scripts\/targetscan7\/Agarwal_;g' | sed -e '86 s;All_cell_lines.AIRs.txt;$ARGV[7];g' | 
sed -e '106 s;1;0;g' |
sed -e '115 s;RNAplfold_in_out;$ARGV[8];g' |  sed -e '118 s;9606;$ARGV[6];g' | sed -e '121 s;qw(10090 10116 13616 8364 9031 9544 9598 9606 9615 9913);("${REF_SPECIES}");g' |
sed -e '631 s;uc($f[3]);$f[3];g' > $2

#!/bin/bash

#    FilTar: Integrating RNA-Seq data to improve microRNA target prediction accuracy in animals
#    Copyright (C) 2019 Thomas Bradley
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>


sed -e '80 s;TA_SPS_by;workflow\/scripts\/targetscan7\/TA_SPS_by;g' $1 | sed -e '83 s;Agarwal_;workflow\/scripts\/targetscan7\/Agarwal_;g' | sed -e '86 s;All_cell_lines.AIRs.txt;$ARGV[7];g' | 
sed -e '106 s;1;0;g' |
sed -e '115 s;RNAplfold_in_out;$ARGV[8];g' |  sed -e '118 s;9606;$ARGV[6];g' | sed -e '121 s;qw(10090 10116 13616 8364 9031 9544 9598 9606 9615 9913);("${REF_SPECIES}");g' |
sed -e '631 s;uc($f[3]);$f[3];g' > $2

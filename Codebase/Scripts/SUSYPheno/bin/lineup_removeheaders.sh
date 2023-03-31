#!/bin/bash
# Usage: ls tables* | lineup_removeheaders.sh [options to lineup.py]
# Ex: ls susy*slha_masses.tlt | lineup_removeheaders.sh -n 2
# The script prints the title line of the first file, and suppresses it for the other files
i=0
(for fn in `cat` ; do i=$[$[i]+1]; if [[ $i == 1 ]]; then head -n 1 ${fn}; fi ; tail -n +2 ${fn} ; done ) | lineup.py $@ 



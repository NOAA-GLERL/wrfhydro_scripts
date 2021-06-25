#!/bin/bash

# print feature_ids and gages to seperate temporary files
ncks --trd -HCs '%i\n' -v link Route_Link.nc > fids.tmp
ncks --trd -HC -v gages Route_Link.nc  | /usr/bin/cut -d\= -f2 | tr -d \' > gages.tmp

# smoooosh together, print only gaged fids and cleanup alignment
paste fids.tmp gages.tmp | awk 'NF>1 {print $0}' | column -t > gages.txt
rm *.tmp




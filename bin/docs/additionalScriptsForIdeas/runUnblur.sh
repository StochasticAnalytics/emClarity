#!/bin/bash

nStart=${1}

nEnd=${2}
fName=${3}

nEnd=$(ls *.mrc | grep -c mrc)
echo $nEnd
for iStack in $(seq 1 $nEnd) ; do


~/unblur_1.0.2/bin/unblur_openmp_7_17_15.exe << eof
frames-${iStack}.mrc
8
frames_${iStack}_avg.mrc
tmp.txt
0.675
no 
no
no
eof

done

#newstack frames_*_avg.mrc gagSP-${fName}-avg.st

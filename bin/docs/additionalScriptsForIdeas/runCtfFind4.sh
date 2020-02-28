#!/bin/bash

ampSpec=${1}
minDef=${2}
maxDef=${3}
exhaustive=${4}

if [ ${exhasutive} != "yes" && ${exhaustive} != "no" ]; then
  echo "specify exhaustive search as yes or no. Try no first."
  exit 1
fi

baseName=$(basename $ampSpec .mrc)
nz=$(header -size ${ampSpec} | awk '{print $3}')

newstack -split 0 -append mrc ${ampSpec} tmpPS_


PIXELSIZE=1.4
KV=300
CS=2.7
AC=0.1
specSize=1024
minRes=24
maxRes=4
defocusStep=30

for i in $(seq 0 $(($nz -1))) ; do 
  j=$(printf "%0.2u" $i)

~/thirdParty/cistem-1.0.0-beta/ctffind --amplitude-spectrum-input << eof &
tmpPS_${j}.mrc
diagnostic_${baseName}_${j}.mrc
$PIXELSIZE
$KV
$CS
$AC
$specSize 
$minRes
$maxRes
$minDef
$maxDef
$defocusStep
no
$exhaustive
no
no
yes
yes
eof
done

wait


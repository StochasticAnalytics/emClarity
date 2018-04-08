#!/bin/bash

baseName=${1}
tltFile=${2}
skipBase=${3}
# assuming these will list in the proper order

newstack ${baseName}??.mrc ${baseName}full.st
rm -f ${baseName}_tmp

ls ${baseName}??.txt | 
  while read a ; do
   echo $a
    tail -n -1 $a |
      awk  '{print (($2-$3)/2)*10^-10,3.141592/180*$4,-1*(($2+$3)/2)*10^-10 }' >> ${baseName}_tmp
  done


awk 'FNR==NR{a[FNR]=$1;b[FNR]=$2;c[FNR]=$3 ;next}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,a[$1],b[$1],$14,c[$1],$16,$17,$18,$19,$20,$21,$22}' ${baseName}_tmp ${tltFile} > ${tltFile}_fixed

#!/bin/bash

# For non-developers git & and merging is overkill. Save a clean copy so that 
# original paths are preserved in custom shell script then push back

emClarityRoot=$(awk '{if(/^emClarity_ROOT=/) print $0}' emClarity)
mcrRoot=$(awk '{if(/^MCR_BASH=/) print $0}' emClarity)

if [[ $mcrRoot && $emClarityRoot ]] ; then
  echo "Found root dirs, moving on."
else
  echo "Did not find root dirs, exiting."
  exit 1
fi

cd ..

mv emClarity emClarity~

git clone --depth=1 https://github.com/bHimes/emClarity.git


awk -v emCR="${emClarityRoot}" '{if(/^emClarity_ROOT=/) print emCR; else print $0}' emClarity/emClarity > emClarity/emClarity.tmp1
awk -v mcrR="${mcrRoot}" '{if(/^MCR_BASH=/) print mcrR; else print $0}' emClarity/emClarity.tmp1 > emClarity/emClarity.tmp2

mv emClarity/emClarity.tmp2 emClarity/emClarity

rm emClarity/emClarity.tmp1

cd emClarity

chmod a=wrx emClarity

emClarity check

echo "Please look at the emClarity check log to make sure all is well."
echo "Then remove the old install which is currently at ../emClarity~"



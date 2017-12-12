#!/bin/bash

# For non-developers git & and merging is overkill. Save a clean copy so that 
# original paths are preserved in custom shell script then push back

emClarityRoot16=$(awk '{if(/^emClarity_ROOT=/) print $0}' emClarity_16b)
mcrRoot16=$(awk '{if(/^MCR_BASH=/) print $0}' emClarity_16b)

emClarityRoot17=$(awk '{if(/^emClarity_ROOT=/) print $0}' emClarity_17b)
mcrRoot17=$(awk '{if(/^MCR_BASH=/) print $0}' emClarity_17b)

if [[ $mcrRoot16 && $emClarityRoot16 ]] ; then
  echo "Found root dirs, moving on."
else
  echo "Did not find root dirs, exiting."
  exit 1
fi

cd ..

echo yes | rm -r emClarity~
mv emClarity emClarity~

git clone --depth=1 https://github.com/bHimes/emClarity.git

# update 16b
awk -v emCR="${emClarityRoot16}" '{if(/^emClarity_ROOT=/) print emCR; else print $0}' emClarity/emClarity_16b > emClarity/emClarity.tmp1
awk -v mcrR="${mcrRoot16}" '{if(/^MCR_BASH=/) print mcrR; else print $0}' emClarity/emClarity.tmp1 > emClarity/emClarity.tmp2

mv emClarity/emClarity.tmp2 emClarity/emClarity_16b

rm emClarity/emClarity.tmp1

# update 17b
awk -v emCR="${emClarityRoot17}" '{if(/^emClarity_ROOT=/) print emCR; else print $0}' emClarity/emClarity_16b > emClarity/emClarity.tmp1
awk -v mcrR="${mcrRoot17}" '{if(/^MCR_BASH=/) print mcrR; else print $0}' emClarity/emClarity.tmp1 > emClarity/emClarity.tmp2

mv emClarity/emClarity.tmp2 emClarity/emClarity_17b

rm emClarity/emClarity.tmp1

cd emClarity

chmod a=wrx emClarity_16b
chmod a=wrx emClarity_17b
#emClarity check

echo "Please look at the emClarity check log to make sure all is well."
echo "Then remove the old install which is currently at ../emClarity~"



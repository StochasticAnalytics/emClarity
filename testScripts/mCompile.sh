#!/bin/bash


shortHead=$(git rev-parse --short HEAD)

mFile=${1}

post="_${shortHead}"

outName="$(basename ${mFile} .m)${post}"

binaryOutName="1_2b_0"
scriptOutName="1_2b_0_19a"



matlab19a -nosplash -nodisplay -nojvm -r "mexCompile ; mcc -m  ${mFile}  -a fitInMap.py  -R -nodisplay -o "$(basename ${mFile} .m)_${binaryOutName}" ; exit" &
      
#I /groups/grigorieff/home/himesb/work/emClarity
    
wait
	rm mccExcludedFiles.log
	rm readme.txt
	rm run_*.sh
	rm requiredMCRProducts.txt

if [ -f emClarity ] ; then
  mv emClarity emClarity~
fi

# Matlab (mcc) complains if ther is an underscore in the name.
#mv emClarity${binaryOutName} emClarity_${binaryOutName}

{

echo '#!/bin/bash'
echo ''
echo '#Please modify this line to point to the text file in your MCR root'
echo '#where you pasted the lines suggested to add to LD_LIBRARY_PATH during install.'
echo '#MCR_BASH="/work/thirdParty/MATLAB/mcr_bash.sh"'
echo 'MCR_BASH=""'
echo ''
echo ''
echo '#Please modify this line to point to the install for emClarity binary'
echo '#emClarity_ROOT=/work/emClarity'
echo 'emClarity_ROOT=""'
echo ''
echo ''
echo ''
echo ''
echo 'if [ -f ${MCR_BASH} ]; then'
echo '  source ${MCR_BASH}'
echo 'else'
echo '  echo "Did not find your mcr_bash file, did you fill it in above?"'
echo '  exit 1'
echo 'fi'
echo ''
echo "if [ ! -f \${emClarity_ROOT}/emClarity_${binaryOutName} ]; then"
echo '  echo "Did not find the binary on the path, did you fill it in above?"'
echo '  exit 1'
echo 'fi'
echo ''
echo "argList="${shortHead} ""
echo 'while [ $# -gt 0 ]; do'
echo '  token=$1'
echo '  argList="${argList} ${token}"'
echo '  shift'
echo 'done'
echo ''
echo "\${emClarity_ROOT}/emClarity_${binaryOutName} \${argList}"
 

} > emClarity_${scriptOutName}

chmod a=wrx emClarity_${scriptOutName}

mkdir -p ../bin
mv emClarity_${scriptOutName} ../bin
#mv emClarity_${binaryOutName} ../../../thirdParty/emClarity
mv emClarity_${binaryOutName} ../bin
cp -rp ../docs ../bin

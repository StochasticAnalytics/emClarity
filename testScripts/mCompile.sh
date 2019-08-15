#!/bin/bash


shortHead=$(git rev-parse --short HEAD)

mFile=${1}

post="_${shortHead}"

outName="$(basename ${mFile} .m)${post}"

major=1
minor=4
bugss=5

#commit=decdc7e
#binaryOutName="testRevert_${commit}"
#scriptOutName="19a_testRevert_${commit}"
binaryOutName="${major}_${minor}_${bugss}"
scriptOutName="mcr_v19a"


matlab19a -nosplash -nodisplay -nojvm -r "mexCompile ; mcc -m  ${mFile}  -a fitInMap.py -a /groups/grigorieff/home/himesb/work/emClarity/mexFiles/compiled/emC_ctfFind -a /groups/grigorieff/home/himesb/work/emClarity/mexFiles/compiled/emC_autoAlign.sh -R -nodisplay -o "$(basename ${mFile} .m)_${binaryOutName}" ; exit" &
      
#I /groups/grigorieff/home/himesb/work/emClarity/mexFiles/compiled/emC_ctffind
    
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
echo 'MCR_BASH=/groups/grigorieff/home/himesb/thirdParty/MTL_MCR_17b/BH_mcr_internal19a.bashrc'
echo ''
echo ''
echo '#Please modify this line to point to the install for emClarity binary'
echo '#emClarity_ROOT=/work/emClarity'
echo 'emClarity_ROOT=/groups/grigorieff/home/himesb/thirdParty/emClarity'
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
mv emClarity_${binaryOutName} ../bin

#cp -rp /nrs/grigorieff/himesb/mayRevertBlank /nrs/grigorieff/himesb/mRevert_${binaryOutName}
#awk -v ON="emClarity_${scriptOutName}" '{if(/^binaryName=/) print "binaryName=/groups/grigorieff/home/himesb/thirdParty/emClarity/"ON; else print $0}' /nrs/grigorieff/himesb/mRevert_${binaryOutName}/runTutorial.sh > /nrs/grigorieff/himesb/mRevert_${binaryOutName}/tmp 
#mv /nrs/grigorieff/himesb/mRevert_${binaryOutName}/tmp /nrs/grigorieff/himesb/mRevert_${binaryOutName}/runTutorial.sh
#chmod a=wrx /nrs/grigorieff/himesb/mRevert_${binaryOutName}/runTutorial.sh 


cp -rp ../docs ../bin
cd ..
zip -r emClarity_${major}.${minor}.${bugss}.zip bin
mv emClarity_${major}.${minor}.${bugss}.zip ~/tmp
#rm -r bin/*






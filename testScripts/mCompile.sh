#!/bin/bash


# NOTE: You will also need to modify your emClarity/mexFiles/mexCompile.m
#   Set the mexPath, and modify the two library linker lines to point at your install of CUDA
#   TODO set up a little configure script to do this and check other deps described below.
# NOTE: You also will need to download the binaries from the emC_dependencies folder on drive
export emC_DEPS="/groups/grigorieff/home/himesb/work/emC_dependencies"

# This grabs the first bit of the commit hash, which then is printed in the logfile
shortHead=$(git rev-parse --short HEAD)

# The is the program you want to compile. Most BH_* can be compiled as standalones, but
# you probably just want the wrapper to then emClarity.m
mFile=${1}

post="_${shortHead}"

outName="$(basename ${mFile} .m)${post}"

# For naming. If you are compiling your own version, use something descriptive in teh
# bugs line. e.g. buggs=5testingFeature
major=1
minor=5
bugs=0
nightly=9

# The final binary, run script and docs folder will be zipped and put in this location
# unless it is NONE then it will be left in the bin dir.
zip_location="${HOME}/tmp"
#zip_location="NONE"


binaryOutName="${major}_${minor}_${bugs}_${nightly}"
scriptOutName="${major}_${minor}_${bugs}_${nightly}_v19a"
#binaryOutName="LTS_fix_${shortHead}"
#scriptOutName="LTS_fix_${shortHead}_v19a"

# You may need to modify this line. 
#     I have "matlab19a" on my path to point to the specific matlab install I want to use.
#     Download the dependencies described in the "statically linked" section here https://github.com/bHimes/emClarity/wiki/Requirements
imodStaticIncludes=""
EMC_ROOT=${HOME}/work/emClarity


#mexCompile ;
matlab19a -nosplash -nodisplay -nojvm -r " mexCompile ; mcc -m  ${mFile} -a fitInMap.py -a ${EMC_ROOT}/alignment/emC_autoAlign -a ${EMC_ROOT}/alignment/emC_findBeads -R -nodisplay -o "$(basename ${mFile} .m)_${binaryOutName}" ; exit" &
      
#I /groups/grigorieff/home/himesb/work/emClarity/mexFiles/compiled/emC_ctffind
    
wait
	rm mccExcludedFiles.log
	rm readme.txt
	rm run_*.sh
	rm requiredMCRProducts.txt

if [ -f emClarity ] ; then
  mv emClarity emClarity~
fi

#Matlab (mcc) complains if ther is an underscore in the name.
#mv emClarity${binaryOutName} emClarity_${binaryOutName}

{

echo '#!/bin/bash'
echo ''
echo '# When this script is invoked, record the PID so that the EMC_tmpDir is deleted'
echo '# even in the event of a crash. (With program script added from EMC_tmpDir.sh)'
echo 'thisPID=$$'
echo ''
echo '#Please modify this line to point to the text file in your MCR root'
echo '#where you pasted the lines suggested to add to LD_LIBRARY_PATH during install.'
echo '#MCR_BASH="/work/thirdParty/MATLAB/mcr_bash.sh"'
echo 'MCR_BASH=/groups/grigorieff/home/himesb/thirdParty/MTL_MCR_17b/BH_mcr_internal19a.bashrc'
echo ''
echo ''
echo '#Please modify this line to point to the install for emClarity binary'
echo '#emClarity_ROOT=/work/emClarity'
echo 'export emClarity_ROOT=/groups/grigorieff/home/himesb/work/emClarity'
echo ''
echo ''
echo '# In order to prevent conflicts with external IMOD or cisTEM binaries, run an '
echo '# your ineractive matlab session through this script.'
echo 'if [[ ${1} == "int" || ${1} == "interactive" ]] ; then'
echo '  echo "Running an interactive matlab session through emClarity"'
echo '  matlabCommand=${2}'
#echo '  resetImodDir=${IMOD_DIR}'
#echo '  unset IMOD_DIR'
#echo '  unset AUTODOC_DIR'
#echo '  export IMOD_DIR=${emClarity_ROOT}/bin/deps'
#echo '  export AUTODOC_DIR=${emClarity_ROOT}/bin/deps/autodoc'
echo '  export IMOD_FORCE_OMP_THREADS=8'
echo '  ${matlabCommand}'
echo '  # Return IMOD to its original state'
#echo '  unset IMOD_DIR'
#echo '  unset AUTODOC_DIR'
#echo '  IMOD_DIR=${resetImodDir}'
echo 'fi'
} > emClarity_${scriptOutName}

cat EMC_tmpDir.sh >> emClarity_${scriptOutName}

{
echo ""
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
#echo 'resetImodDir=${IMOD_DIR}'
#echo 'unset IMOD_DIR'
#echo 'export IMOD_DIR=${emClarity_ROOT}/bin/deps'
echo "\${emClarity_ROOT}/emClarity_${binaryOutName} \${argList}"
#echo 'export IMOD_DIR=${resetImodDir}'
 

} >> emClarity_${scriptOutName}

chmod a=wrx emClarity_${scriptOutName}

# Collect the emClarity binary dependencies
mkdir -p ../bin
mkdir -p ../bin/deps
#mkdir -p ../bin/deps/autodoc
#mkdir -p ../bin/deps/com
#mkdir -p ../bin/deps/bin
#mkdir -p ../bin/deps/bin/realbin
#cp -ru ${emC_DEPS}/deps/bin/realbin/* ../bin/deps/bin/realbin
#cp -ru ${emC_DEPS}/deps/com/* ../bin/deps/com

#cp -ru ${emC_DEPS}/deps/VERSION ../bin/deps
#cp -ru ${emC_DEPS}/deps/imodDeps.txt ../bin/deps
cp -ru ${emC_DEPS}/deps/cisTEMDeps.txt ../bin/deps
#cd ../bin/deps/autodoc
#cat ../imodDeps.txt | while read dep ; do
#  cp -u ${emC_DEPS}/deps/autodoc/${dep}.adoc .
#  ln -sf ${dep}.adoc emC_${dep}.adoc
#done
cd ${EMC_ROOT}/testScripts
cat ../bin/deps/cisTEMDeps.txt | while read dep ; do
  cp -u ${emC_DEPS}/emC_${dep} ../bin/deps
done

mv emClarity_${scriptOutName} ../bin
mv emClarity_${binaryOutName} ../bin

#cp -rp /nrs/grigorieff/himesb/mayRevertBlank /nrs/grigorieff/himesb/mRevert_${binaryOutName}
#awk -v ON="emClarity_${scriptOutName}" '{if(/^binaryName=/) print "binaryName=/groups/grigorieff/home/himesb/thirdParty/emClarity/"ON; else print $0}' /nrs/grigorieff/himesb/mRevert_${binaryOutName}/runTutorial.sh > /nrs/grigorieff/himesb/mRevert_${binaryOutName}/tmp 
#mv /nrs/grigorieff/himesb/mRevert_${binaryOutName}/tmp /nrs/grigorieff/himesb/mRevert_${binaryOutName}/runTutorial.sh
#chmod a=wrx /nrs/grigorieff/himesb/mRevert_${binaryOutName}/runTutorial.sh 


cp -rp ../docs ../bin
cd ..
if [[ ${zip_location} != "NONE" ]]; then
  zip -r emClarity_${major}.${minor}.${bugs}.${nightly}.zip bin
  mv emClarity_${major}.${minor}.${bugs}.${nightly}.zip ${zip_location}
fi






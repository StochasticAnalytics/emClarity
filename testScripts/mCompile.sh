#!/bin/bash


# NOTE: You will also need to modify your emClarity/mexFiles/mexCompile.m
#   Set the mexPath, and modify the two library linker lines to point at your install of CUDA
#   TODO set up a little configure script to do this and check other deps described below.
# NOTE: You also will need to download the binaries from the emC_dependencies folder on drive
export emC_DEPS="/sa_shared/software/emClarity_1.6.1.0/bin/deps"
EMC_ROOT=/sa_shared/software/ # This is just for convenience, when you build for yourself
############ lines 4-5 in mexFiles/mexCompile
#mexPATH = '/groups/grigorieff/home/himesb/work/emClarity/mexFiles/';
#CUDA_LIB = '-L/groups/grigorieff/home/himesb/thirdParty/cuda-10.0/lib64'   ... % NOTE if you leave a space at the end of this string, MATLAB does not parse the option correctly (which wouldn't matter in a normal compile line!)
############
### line 7 in testScripts/emClarity.m
#compiled_PATH='/groups/grigorieff/home/himesb/work/emClarity';
###################

# This is the version of matlab you will end up compiling with.
MATLAB_FOR_COMIPLING=matlab

# This grabs the first bit of the commit hash, which then is printed in the logfile
shortHead=$(git rev-parse --short HEAD)

# The is the program you want to compile. Most BH_* can be compiled as standalones, but
# you probably just want the wrapper to then emClarity.m
mFile=${1}

post="_${shortHead}"

# First make sure the file exists
if [[ ! -f ${mFile} ]]; then
  echo "Could not find ${mFile}"
  exit 1
fi
# check the extension
if [[ $mFile != *.m ]]; then
  echo "Please provide a matlab file to compile"
  exit 1
fi

outName="$(basename ${mFile} .m)${post}"

# For naming. If you are compiling your own version, use something descriptive in teh
# bugs line. e.g. buggs=5testingFeature
major=1
minor=7
bugs=0
nightly=7

EMC_ROOT=${EMC_ROOT}/emClarity_${major}.${minor}.${bugs}.${nightly}

# The final binary, run script and docs folder will be zipped and put in this location
# unless it is NONE then it will be left in the bin dir.
zip_location="${HOME}/tmp"
#zip_location="NONE"


binaryOutName="${major}_${minor}_${bugs}_${nightly}"
scriptOutName="${major}_${minor}_${bugs}_${nightly}_v23a"
#binaryOutName="LTS_fix_${shortHead}"
#scriptOutName="LTS_fix_${shortHead}_v19a"

# You may need to modify this line. 
#     I have "matlab19a" on my path to point to the specific matlab install I want to use.
#     Download the dependencies described in the "statically linked" section here https://github.com/bHimes/emClarity/wiki/Requirements
imodStaticIncludes=""




${MATLAB_FOR_COMIPLING} -nosplash -nodisplay -nojvm -r " mexCompile ; mcc -m  ${mFile} -a fitInMap.py -a ../alignment/emC_autoAlign -a ../alignment/emC_findBeads -a ../metaData/BH_checkInstall -R -nodisplay -o "$(basename ${mFile} .m)_${binaryOutName}" ; exit" &
          
wait

rm mccExcludedFiles.log
rm readme.txt
rm run_*.sh
rm requiredMCRProducts.txt
rm unresolvedSymbols.txt
rm includedSupportPackages.txt

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
echo ''
echo '# Note you no longer need to modify this line inside the singularity container:'
echo 'MCR_BASH=/sa_shared/software/matlab2023/MATLAB/R2023a/runtime/glnxa64:/sa_shared/software/matlab2023/MATLAB/R2023a/bin/glnxa64:/sa_shared/software/matlab2023/MATLAB/R2023a/sys/os/glnxa64'
echo ''
echo ''
echo '#Please modify this line to point to the install for emClarity binary'
echo '#emClarity_ROOT=${HOME}/emC_builds'
echo "export emClarity_ROOT=${EMC_ROOT}"
echo 'export LD_LIBRARY_PATH=${emClarity_ROOT}/lib:${MCR_BASH}:${LD_LIBRARY_PATH}'
echo ''

} > emClarity_${scriptOutName}

cat EMC_tmpDir.sh >> emClarity_${scriptOutName}

{
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

} >> emClarity_${scriptOutName}

{
echo ""
echo ''
echo "if [ ! -f \${emClarity_ROOT}/bin/emClarity_${binaryOutName} ]; then"
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
echo "\${emClarity_ROOT}/bin/emClarity_${binaryOutName} \${argList}"
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

rm -rf ../emClarity_${major}.${minor}.${bugs}.${nightly}
mkdir ../emClarity_${major}.${minor}.${bugs}.${nightly}
cp -rp ../docs ../emClarity_${major}.${minor}.${bugs}.${nightly}
cp -rp ../lib ../emClarity_${major}.${minor}.${bugs}.${nightly}

cd ../emClarity_${major}.${minor}.${bugs}.${nightly}
mkdir bin
cp -rp ../bin/deps ./bin
cd bin
cp ${EMC_ROOT}/testScripts/.bashrc .
cp ../../bin/emClarity_${scriptOutName} emClarity_${scriptOutName}
cp ../../bin/emClarity_${binaryOutName} emClarity_${binaryOutName}
ln -s emClarity_${scriptOutName} emClarity
cd ../../

if [[ ${zip_location} != "NONE" ]]; then
  zip -r --symlinks emClarity_${major}.${minor}.${bugs}.${nightly}.zip ./emClarity_${major}.${minor}.${bugs}.${nightly}
  mv emClarity_${major}.${minor}.${bugs}.${nightly}.zip ${zip_location}
fi

rm -rf emClarity_${major}.${minor}.${bugs}.${nightly}





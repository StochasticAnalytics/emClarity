#!/bin/bash


# NOTE: You will also need to modify your emClarity/mexFiles/mexCompile.m
#   Set the mexPath, and modify the two library linker lines to point at your install of CUDA
#   TODO set up a little configure script to do this and check other deps described below.
# NOTE: You also will need to download the binaries from the emC_dependencies folder on drive
export emC_DEPS="/sa_shared/software/emClarity_1.6.1.0/bin/deps"
EMC_COMPILED_ROOT=/sa_shared/software/ # This is just for convenience, when you build for yourself
############ lines 4-5 in mexFiles/mexCompile
#mexPATH = '/groups/grigorieff/home/himesb/work/emClarity/mexFiles/';
#CUDA_LIB = '-L/groups/grigorieff/home/himesb/thirdParty/cuda-10.0/lib64'   ... % NOTE if you leave a space at the end of this string, MATLAB does not parse the option correctly (which wouldn't matter in a normal compile line!)
############
### line 7 in testScripts/emClarity.m
#compiled_PATH='/groups/grigorieff/home/himesb/work/emClarity';
###################

# This is the version of matlab you will end up compiling with.
MATLAB_FOR_COMPILING=matlab

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
nightly=15

binaryOutName="${major}_${minor}_${bugs}_${nightly}"
scriptOutName="${major}_${minor}_${bugs}_${nightly}_v23a"

EMC_VERSION=emClarity_${major}.${minor}.${bugs}.${nightly}

EMC_COMPILED_DIRNAME=${EMC_COMPILED_ROOT}/${EMC_VERSION}

# The final binary, run script and docs folder will be zipped and put in this location
# unless it is NONE then it will be left in the bin dir.
zip_location="${HOME}/tmp"
#zip_location="NONE"



#binaryOutName="LTS_fix_${shortHead}"
#scriptOutName="LTS_fix_${shortHead}_v19a"

# You may need to modify this line. 
#     I have "matlab19a" on my path to point to the specific matlab install I want to use.
#     Download the dependencies described in the "statically linked" section here https://github.com/bHimes/emClarity/wiki/Requirements
imodStaticIncludes=""

# Generally only enabled to speed up debugging of compilation steps unrelated to the cudaMex files which are otherwise always compiled first.
skipMex=0
if [[ $skipMex -eq 0 ]]; then
  mexCompile="mexCompile ;"
else
  mexCompile=""
fi

# NOTE: warnings are disabled to ensure that failed builds are caught. Ideally they would be addressed and removed.
${MATLAB_FOR_COMPILING} -nosplash -nodisplay -nojvm -r " ${mexCompile} mcc -w disable -w off -m  ${mFile} -a fitInMap.py -a ../alignment/emC_autoAlign -a ../alignment/emC_findBeads -a ../metaData/BH_checkInstall -R -nodisplay -o "$(basename ${mFile} .m)_${binaryOutName}" ; exit" &
          
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
echo "export emClarity_ROOT=${EMC_COMPILED_DIRNAME}"
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


rm -rf ../bin/${EMC_VERSION}
mkdir ../bin/${EMC_VERSION} ../bin/${EMC_VERSION}/bin ../bin/${EMC_VERSION}/bin/deps

mv emClarity_${scriptOutName} ../bin/${EMC_VERSION}/bin
mv emClarity_${binaryOutName} ../bin/${EMC_VERSION}/bin

cp -ru ${emC_DEPS}/deps/cisTEMDeps.txt ../bin/${EMC_VERSION}/bin/deps

cat ../bin/${EMC_VERSION}/bin/deps/cisTEMDeps.txt | while read dep ; do
  echo "Copying ${dep} to ../bin/${EMC_VERSION}/bin/deps"
  cp -u ${emC_DEPS}/emC_${dep} ../bin/${EMC_VERSION}/bin/deps
done

cp -rp ../docs ../bin/${EMC_VERSION}
cp -rp ../lib ../bin/${EMC_VERSION}

cp .bashrc ../bin/${EMC_VERSION}
cd ../bin/${EMC_VERSION}/bin

ln -s emClarity_${scriptOutName} emClarity
cd ../../

if [[ ${zip_location} != "NONE" ]]; then
  zip -r --symlinks ${EMC_VERSION}.zip ./${EMC_VERSION}
  mv ${EMC_VERSION}.zip ${zip_location}
fi

rm -rf ${EMC_VERSION}





#!/bin/bash


BH_wrkdir=${PWD}

if [[ ${1} -eq 0 ]] ; then
  shift
  depend=""
else
  depend="-W depend=afterany:${1}"
  shift
fi

# convert the argument list to a string.
	programName=${1}
	shift
  paramName=${1}
	# Second arg is param file, save a copy.
  argList=
  while [ $# -gt 0 ]; do
      token=$1
      argList="${argList} ${token}" 
      shift
  done
echo $argList
mkdir -p paramRecord

dT="$(echo $(date) | awk '{print $4}' | awk -F ":" '{print $1"_"$2"_"$3}')"
jobName=job_${programName}_${dT}
{
	echo '#!/bin/bash'
	echo ''
	#echo "#PBS -N job_${programName}"
	echo '#PBS -l nodes=1:ppn=28:gpus=4'
	echo '#PBS -l walltime=04:00:00:00'
	echo '#PBS -V'
	echo '#PBS -j oe'
	echo '#PBS -m e'
#	echo '#PBS -M yourEmailAT@wherever.com'
	echo ''
	echo 'source ${HOME}/MATLAB_BH_v15a/mcr2015a_BH.bash'
        echo 'module load cuda75/toolkit/7.5.18'
	echo 'module load cryoem/IMOD/4.8.57'
	echo ''
  
	echo "${programName} ${argList}"
  echo "logFile=\$(ls -rth *.log | tail -n -1)"
  echo "mv \${logFile} ${dT}_\${logFile}"
 # echo "kms=\$(ls -rth cycle001*Kms.mrc | tail -n -1)"
 # echo "mv \${kms} ${paramName}.mrc"
  echo "rm BH*.o*"
  echo "rm ${jobName}"
	echo ''
}	> ${jobName}
	

echo "qsub ${depend} -d ${BH_wrkdir}  job_${programName}" 
sent=$( qsub ${depend} -d ${BH_wrkdir} ${jobName} -N BH_${dT} | awk -F "." '{ print $1 }' )
echo $sent
echo $sent > lastSubmissionIDX.txt


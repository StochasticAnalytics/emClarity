# This is added in the middle of the emClarity run script. It does the following.
# 1 - check for temporary file systems so that some disk i/o can instead be kept in memory.
# 2 - minimally the MCR_CACHE_ROOT is set there (formerly .mcrCarch_9.3 which is default)
#     This avoids the problem where parallel pools are dying due to either slow filesystems
#     or conflicts with other parpools from competing nodes.
# 3 - Fall back to $TMPDIR or /tmp if no tmpfs is found. 
# 4 - Some emClarity funcitons will use this space too iff there is enough memory available.
#     Note that the mem avail is queried at runtime, so this could be a problem TODO

################################################################################
################################################################################
#### NEW section to automatically create and delete a tmp cache dir.

# Default to /tmp if the following doesn't workspace

MAX_MEM=0
MAX_FS=''
EMC_CACHE_DIR=''
EMC_CACHE_MEM=''

while read target; do

  if [ -w $target ] ; then 

    this_mem=$(df --output=avail ${target} | tail -n -1)
    if [[ $this_mem && $this_mem -gt $MAX_MEM ]] ; then 
      MAX_MEM=${this_mem}
      MAX_FS=${target}
    fi

  fi

done <<< "$(findmnt -t tmpfs --output=TARGET )"

if [ $MAX_FS ] ; then
  echo "Found tmpfs $MAX_FS with max mem $MAX_MEM bytes"
  echo "using this for MCR_CACHE_ROOT"
  export EMC_CACHE_DIR=$MAX_FS
else
  if [ $TMPDIR ] ; then 
    export EMC_CACHE_DIR=${TMPDIR}
  else
    export EMC_CACHE_DIR="/tmp"
  fi
  echo "Did not find a suitable tmpfs, using default $MAX_FS"
fi



TMP_DIR_ID=emC_tmp_${RANDOM}

# Set the MCR cache their
unset MCR_CACHE_ROOT
export MCR_CACHE_ROOT="${EMC_CACHE_DIR}/${TMP_DIR_ID}"

mkdir -p ${MCR_CACHE_ROOT}

# Max time for a job in minutes. Default 7 days (10080)
MAX_TIME=10080
# Startup a background process that will make sure the tmp dir is removed
# even on crash
{
  echo "elapsed_time=0"
  echo "while [[ \$(ps -q ${thisPID} -o comm=) ]] ; do"
  echo "  sleep 60; if [[ \${elapsed_time} > $MAX_TIME ]] ; then break; fi"
  echo "  elapsed_time=$((${elapsed_time} + 1))"
  echo "done"
  echo "rm -rf ${MCR_CACHE_ROOT}"

} > ${MCR_CACHE_ROOT}/emC_memClean.sh
chmod a=wrx ${MCR_CACHE_ROOT}/emC_memClean.sh
${MCR_CACHE_ROOT}/emC_memClean.sh &

export EMC_CACHE_MEM=$(echo $MAX_MEM | awk '{print $0/1024/1024}')
echo -e "\n##############\n"
echo -e "\tCreated a tmp MCR_CACHE locally at ${MCR_CACHE_ROOT}"
echo -e "\tAvailable mem is $EMC_CACHE_MEM GB"
echo -e "\tStarted emC_memClean to make sure it removed in the event of a crash."
echo -e "\tDefault kill time is a max of $MAX_TIME minutes"
echo -e "\n##############\n"
################################################################################
################################################################################

echo "Checking the installation of emClarity and some system variables." > emClarity_checkInstall.txt

{
echo -e "\nRun on : $data" 
echo "Recording variables that SHOULD be set in the users profile prior to entering emClarity" 
} >>  emClarity_checkInstall.txt

# These should be set in the users shell
for i in CUDA_VISIBLE_DEVICES CUDA_BIN_PATH CUDA_HOME CUDA_CACHE_MAXSIZE CUDA_CACHE_DISABLE CUDA_CACHE_PATH MATLAB_SHELL CUDA_VERSION=$(nvcc --version | awk -F "release" /release/'{print $2}') IMOD_DIR ; do

  echo "${i}:$(echo $i)" >>  emClarity_checkInstall.txt

done




SHELL: This describes the shell that will be interpreting any commands you type in. In most cases, this will be bash by default, but other values can be set if you prefer other options.
TERM: This specifies the type of terminal to emulate when running the shell. Different hardware terminals can be emulated for different operating requirements. You usually won’t need to worry about this though.
USER: The current logged in user.
PWD: The current working directory.
PATH: A list of directories that the system will check when looking for commands. When a user types in a command, the system will check directories in this order for the executable.
LANG: The current language and localization settings, including character encoding.
HOME: The current user’s home directory.


# These are all set in the emClarity run script
EMC_CACHE_DIR
EMC_CACHE_MEM
EMC_ROOT
MCR_CACHE_ROOT
emClarity_ROOT
IMOD_FORCE_OMP_THREADS
emC_DEPS

# These are set in the emClarity binary
MATLAB_SHELL
EMC_AUTOALIGN
EMC_FINDBEADS

emC_cisTEMDepPath=sprintf('%s/bin/deps',emClarity_ROOT)
emC_cisTEM_deps = importdata(sprintf('%s/cisTEMDeps.txt',emC_cisTEMDepPath));
for iDep = 1:length(emC_cisTEM_deps)
  setenv(sprintf('EMC_%s',upper(emC_cisTEM_deps{iDep})), sprintf('%s/emC_%s',emC_cisTEMDepPath,emC_cisTEM_deps{iDep}));
end
  setenv('CUDA_VISIBLE_DEVICES',cudaStart);


echo "Checking the installation of emClarity and some system variables." > emClarity_checkInstall.txt

{
echo -e "\nRun on : $data" 
echo "Recording variables that SHOULD be set in the users profile prior to entering emClarity" 
} >>  emClarity_checkInstall.txt

echo "CUDA_VISIBLE_DEVICES:$CUDA_VISIBLE_DEVICES" >>  emClarity_checkInstall.txt
echo "CUDA_HOME:$(echo $CUDA_HOME)" >>  emClarity_checkInstall.txt
echo "CUDA_BIN_PATH:$(echo $CUDA_BIN_PATH)" >>  emClarity_checkInstall.txt
echo "CUDA_CACHE_MAXSIZE:$(echo $CUDA_CACHE_MAXSIZE)" >>  emClarity_checkInstall.txt
echo "CUDA_CACHE_DISABLE:$(echo $CUDA_CACHE_DISABLE)" >>  emClarity_checkInstall.txt
echo "CUDA_CACHE_PATH:$(echo $CUDA_CACHE_PATH)" >>  emClarity_checkInstall.txt
echo "MATLAB_SHELL:$(echo $MATLAB_SHELL)" >>  emClarity_checkInstall.txt
echo "IMOD_DIR:$(echo $IMOD_DIR)" >>  emClarity_checkInstall.txt


echo -e "\n\nRecording some system variables that might be helpful in troubleshooting"  >>  emClarity_checkInstall.txt


echo "SHELL:$(echo $SHELL)" >>  emClarity_checkInstall.txt
echo "TERM:$(echo $TERM)" >>  emClarity_checkInstall.txt
echo "USER:$(echo $USER)" >>  emClarity_checkInstall.txt
echo "PWD:$(echo $PWD)" >>  emClarity_checkInstall.txt
echo "PATH:$(echo $PATH)" >>  emClarity_checkInstall.txt
echo "LANG:$(echo $LANG)" >>  emClarity_checkInstall.txt
echo "HOME:$(echo $HOME)" >>  emClarity_checkInstall.txt



echo -e "\n\nRecording some variables that are set in the emClarity run script"  >>  emClarity_checkInstall.txt


echo "EMC_CACHE_DIR:$(echo $EMC_CACHE_DIR)" >>  emClarity_checkInstall.txt
echo "EMC_CACHE_MEM:$(echo $EMC_CACHE_MEM)" >>  emClarity_checkInstall.txt
echo "MCR_CACHE_ROOT:$(echo $MCR_CACHE_ROOT)" >>  emClarity_checkInstall.txt
echo "emClarity_ROOT:$(echo $emClarity_ROOT)" >>  emClarity_checkInstall.txt
echo "IMOD_FORCE_OMP_THREADS:$(echo $IMOD_FORCE_OMP_THREADS)" >>  emClarity_checkInstall.txt



echo -e "\n\nRecording some variables that are set in the emClarity binary"  >>  emClarity_checkInstall.txt

echo "MATLAB_SHELL:$(echo $MATLAB_SHELL)" >>  emClarity_checkInstall.txt
echo "EMC_AUTOALIGN:$(echo $EMC_AUTOALIGN)" >>  emClarity_checkInstall.txt
echo "EMC_FINDBEADS:$(echo $EMC_FINDBEADS)" >>  emClarity_checkInstall.txt



#!/bin/bash


baseName=${1}
modFile="bin10/${1}_bin10.mod"
modBin=10 # this could be changed, but works well in most cases and is ill advised.
modThick=300 # this could be changed, particularly if your sample is thicker
             # than 2000 pixels. But then it is probably a bad candidate for
             # high resolution tomography anyhow ( assuming your pixel size > 1Ang) 

wrkDir=$(pwd)

if [[  ! -f aliStacks/${baseName}_ali1.fixed ]] ; then
  echo -e "\n\nDid not find the aligned, fixed stack at ${baseName}_ali1.fixed\n\n"
  exit 1
fi

# Get the size of the stack
sizeSt=($(header -size aliStacks/${baseName}_ali1.fixed))

#optional
if [[ -f fixedStacks/${baseName}.local ]] ; then
  echo -e "Using local alignments, and trusting the header is set properly.\n"
  echo -e "imod will compare column 8 and the x,y,z pixel sizes to scale local.\n"
  echo ""
  head -n 1 fixedStacks/${baseName}.local
  echo ""
  echo $(header -pixel aliStacks/${baseName}_ali1.fixed)
  localFILE="-LOCALFILE ${wrkDir}/fixedStacks/${baseName}.local \\"
else
  localFILE=""
  echo "WARNING - no local files found at fixedStacks/${baseName}.local"
  echo -e "\nproceeding anyway.\n"
fi


model2point -float -contour ${modFile} mod.tmp
sort -k 1 -g mod.tmp > ${modFile}.txt
rm mod.tmp

numParticles=($(wc ${modFile}.txt | tail -n -1 | awk '{print $1/6,$1%6}'))
echo ${numParticles[1]}
echo ${numParticles[0]}
if [[ ${numParticles[1]} -gt 0 ]] ; then
  echo -e "\n\nThe number of contours (${numParticles[0]}) is not a multiple of 6\n\n"
  exit 1
fi



mkdir -p recon

# Write the basename and num particles once followed by the recon parameters

 echo "${baseName}" > recon/${baseName}_recon.coords
 echo "${numParticles[0]}" >> recon/${baseName}_recon.coords
 

# For gag particles, pick right at the apparent edge, and use a pad factor of
# 1.2 - For other samples where you are just reducing the area, use a pad factor of
# 1.0 or at least check that no out of bounds are included (add auto check someday)
padFact=1.07
for iPart in $(seq 1 ${numParticles[0]}) ; do 
  c=[]  

  unset xLt xHt yLt yHt zLt zHt WIDTH Height yLow yTop oX THICKNESS oZ
  
  pN=$(echo "print (($iPart-1)*6)+1" | python)
  xLt=$(awk -v T=$(($pN + 0)) -v B=${modBin} '{if(FNR==T) print $2*B }' ${modFile}.txt)
  xHt=$(awk -v T=$(($pN + 1)) -v B=${modBin} '{if(FNR==T) print $2*B  }' ${modFile}.txt)
  yLt=$(awk -v T=$(($pN + 2)) -v B=${modBin} '{if(FNR==T) print $3*B }' ${modFile}.txt)
  yHt=$(awk -v T=$(($pN + 3)) -v B=${modBin} '{if(FNR==T) print $3*B }' ${modFile}.txt)
  zLt=$(awk -v T=$(($pN + 4)) -v B=${modBin} '{if(FNR==T) print $4*B  }' ${modFile}.txt) 
  zHt=$(awk -v T=$(($pN + 5)) -v B=${modBin} '{if(FNR==T) print $4*B  }' ${modFile}.txt)

  # Explicitly check min and max
  xL=$(echo "print min($xLt,$xHt)" | python)
  xH=$(echo "print max($xLt,$xHt)" | python)
  yL=$(echo "print min($yLt,$yHt)" | python)
  yH=$(echo "print max($yLt,$yHt)" | python)
  zL=$(echo "print min($zLt,$zHt)" | python)
  zH=$(echo "print max($zLt,$zHt)" | python)


  echo "xl $xL xh $xH yl $yL yh $yH zl $zL zh $zH"
  WIDTH=$(echo "print int((${padFact}*($xH-$xL)))" | python)
  Height=$(echo "print int((${padFact}*($yH-$yL)))" | python)
  echo "$WIDTH $Height"
  yLow=$(echo "print int($yL)" | python)
  yTop=$(echo "print int(($yLow + $Height - 1 ))" | python)
  echo "$yLow $yTop"
  oX=$(echo "print int((${sizeSt[0]}/2.0 - (${xL}+(${WIDTH}/2.0))))" | python)
  echo "$oX"
  THICKNESS=$(echo "print int((${padFact}*(($zH - $zL)**2)**0.5))" | python)
  echo "$THICKNESS"
  oZ=$(echo "print -1*int((($modThick*$modBin)/2.0- ($zL+$THICKNESS/2.0)))" | python)
  echo "$oZ"
  oY=$(echo "print int(0.5*(${yLow} + ${yTop} - 1) - ${sizeSt[1]}/2.0 ) " | python)
  
# THICKNESS=$(($THICKNESS+140))

# These files will be pulled into the metadata and used if you decide to 
# do particle polishing, i.e. refine the tilt series geometry using 
# subTomogram based Fiducial markers. After each round the ali1 will change
# accordingly ali2 ...etc
echo "$WIDTH $Height $THICKNESS $((-1*$oX)) $oY $oZ" > recon/${baseName}_${iPart}_recon.txt
{
 
  echo "${WIDTH}"
  echo "$yLow"
  echo "$yTop"
  echo "$THICKNESS"
  echo "$oX"
  echo "$oZ"


} >> recon/${baseName}_recon.coords

# For the initial cycle, prior to template matching and any subTomo work
# Equivalent bash scripts are written in the same directory. Call these directly 
# And if you decide any tomos aren't worth keeping, the .coords will be ignored
# if you delete the .path file that is in the template matching results directory
# called convmap.
{

  echo '#!/bin/bash'
  echo ''
  echo "tilt -input ${wrkDir}/ctfStacks/${baseName}_ali1_ctf.fixed \\"
  echo "     -TILTFILE ${wrkDir}/fixedStacks/${baseName}.tlt \\"
  echo "     -output ${wrkDir}/recon/${baseName}_${iPart}.rec \\"
  echo "     -RADIAL 0.475,0.05 \\"
  echo "     -UseGPU 0 \\"
  echo "     -WIDTH $WIDTH \\"
  echo "     -SLICE $yLow,$yTop \\"
  echo "     -THICKNESS $THICKNESS \\"
  echo "     -RotateBy90 \\"
  if [[ ${localFILE} ]] ; then
  echo "     ${localFILE}"
  fi
  echo "     -SHIFT $oX,$oZ"



} > recon/${baseName}_${iPart}_recon.sh

chmod a=wrx recon/${baseName}_${iPart}_recon.sh
pN=$(($pN+5))
 
done

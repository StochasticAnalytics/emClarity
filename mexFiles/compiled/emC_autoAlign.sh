#!/bin/bash

#Some fixed values that may be edited, but seem to be generally good:

  # Resolution cutoff is best to be at or slightly lower resolution than the first zero of the CTF
  # TODO add input option for defocus value estimate, and calculate the first zero automatically.
  # This would also facilitate a quick handedness check prior to full ctf estimation.
  RESOLUTION_CUTOFF=15
  LOW_RES_CUTOFF=1800
  # Angstrom, edge width
  PATCH_SIZE=1200 
  # The run-time vs improved accuracy is much worse for small pixel sizes. (ang/pix)
  MAX_SAMPLING_RATE=3.0
  MIN_SAMPLING_RATE=10.0
  ITERATIONS_PER_BIN=3
  TILT_OPTION=5
  MAG_OPTION=5
  tiltAngleOffset=0.0


  # It really seems like the high-pass filter in tiltxcorr is not being used
  FORCE_HIGH_PASS=false

# TODO make the default run etomo on first iteration for inspection. Add flag to override for clean datasets.

# Check for command line input, if not, as the user for a file and pixel size.
# TODO add a wild card option or directory check to run automatically for all stacks
foundCmdLineArgs=false



inp=${1}
pixelSize=${2}
tiltAxisRotation=${3}
binHigh=${4}
binLow=${5}
binInc=${6}
nX=${7}
nY=${8}
nZ=${9}
ext=${10}

echo $pwd

iter=1;
echo " $binHigh $binInc $binLow"
echo "$(seq $binHigh $binInc $binLow)"
for iBin in $(seq $binHigh $binInc $binLow) ; do 
#if [[ $iter -gt 3 ]] ; then exit 1; fi

  lpCut=$(echo "print(${iBin}*${pixelSize}/${RESOLUTION_CUTOFF})" | python)
  hpCut=$(echo "print(${iBin}*${pixelSize}/${LOW_RES_CUTOFF})" | python)
  echo -e "\n\n lp cutoff $RESOLUTION_CUTOFF , $lpCut, hp cutoff $hpCut\n\n"

  for iBinIter in $(seq 1 $ITERATIONS_PER_BIN) ; do

    echo " ${iBin}"
    if [[ ${iter} -lt 4 ]] ; then

      ptSize=$(echo "print(int(${nX}/($iter*${iBin})))" | python)
      ptBorder=64
      ptOverlap=0.5
      # For each binning the pixel shifts are limited to 50 A * iBin
#    elif [[ ${iBin} -eq 2 ]]  ; then
#      ptSize=$(echo "print int(${PATCH_SIZE}/(${iBin}*${pixelSize}))" | python)
#      ptBorder=128
#      ptOverlap=0.5
#      limitShifts=$(echo "print int(150.0/($iter*${pixelSize}))"  | python)
    else
      ptSize=$(echo "print(int(${nX}/($iBin*3.)))" | python)
      ptBorder=64
      ptOverlap=0.5

    fi

    firstIterShiftLimit=10
    if [[ $iter -ge 2 ]] ; then
      limitShifts=$(echo "print(int(${firstIterShiftLimit}/$iter)+1)" | python)
    else
      limitShifts=${firstIterShiftLimit}
    fi
 
    wDir="round_${iter}"
    pName="INP-${iter}"

    mkdir -p ${wDir}


    cd ${wDir}


    if [[ ${iter} -eq 1 ]] ; then
      # open up etomo to run through to visually inspect ccd removal and get the 
      # alignment to be sane. assuming tilts in header or manually extracted.
      ln -s ../../${inp}${ext} ${pName}.st
      ln -s ../../${inp}.rawtlt ${pName}.rawtlt
#      if [[ -f ../../$(basename ${inp} .mrc).rawtlt ]] ; then
#        inp=$(basename ${inp} .mrc)
#        ln -s ../../${inp}.rawtlt ${pName}.rawtlt        
#      elif [[ ../../$(basename ${inp} .st).rawtlt ]] ; then
#        inp=$(basename ${inp} .st)
#        ln -s ../../${inp}.rawtlt ${pName}.rawtlt    
#      elif [[ ../../$(basename ${inp} .fixed).rawtlt ]] ; then
#        inp=$(basename ${inp} .fixed)
#        ln -s ../../${inp}.rawtlt ${pName}.rawtlt 
#      else 
#        echo "Could not find a *.rawtlt file. Is it there? is your stack a .mrc or .st?"
#        exit 1
#      fi

      
#      ccderaser \
#        -InputFile	${pName}.st \
#        -OutputFile	../../fixedStacks/${inp}.fixed \
#        -FindPeaks	\
#        -PeakCriterion	7. \
#        -DiffCriterion	8. \
#        -BigDiffCriterion	14. \
#        -GiantCriterion	12. \
#        -ExtraLargeRadius	8. \
#        -GrowCriterion	4. \
#        -EdgeExclusionWidth	4 \
#        -PointModel	${pName}_peak.mod \
#        -MaximumRadius	2.1 \
#        -AnnulusWidth	2.0 \
#        -XYScanSize	100 \
#        -ScanCriterion	3. \
#        -BorderSize	2 \
#        -PolynomialOrder	2

#      newstack -xf ../../preRotXf.xf ../../fixedStacks/${inp}.fixed ../../fixedStacks/${inp}.fixed.rot

     
#      if [[ $FORCE_HIGH_PASS == true ]] ; then
#        echo -e "\n\nForcing high-pass filter\n\n"
#        sleep 2
#        maxCut=$(echo "print(${pixelSize}/${LOW_RES_CUTOFF})" | python)
#        mv ../../fixedStacks/${inp}.fixed.rot ../../fixedStacks/${inp}.fixedNP.rot
#        clip filter -h ${maxCut} ../../fixedStacks/${inp}.fixedNP.rot ../../fixedStacks/${inp}.fixed.rot
#      fi

      tiltxcorr \
        -InputFile	../../fixedStacks/${inp}.fixed.rot \
        -OutputFile	${pName}.prexf \
        -TiltFile	${pName}.rawtlt \
        -RotationAngle	0.0 \
        -FilterSigma1	0.001 \
        -FilterRadius2	0.5 \
        -FilterSigma2	0.05 \
        -CumulativeCorrelation \
        -AngleOffset ${tiltAngleOffset}
        	


      xftoxg -nfit 0 ${pName}.prexf ${pName}.inpXF

      cp ${pName}.rawtlt ${pName}.tlt
      
      cd ..
      iter=$(($iter+1))

      continue
    else

      # TODO shouldn't this be set in the first iter? Make it persist.
      if [[ -f ../../$(basename ${inp} .mrc).rawtlt ]] ; then
        inp=$(basename ${inp} .mrc)
      elif [[ ../../$(basename ${inp} .st).rawtlt ]] ; then
        inp=$(basename ${inp} .st)
      fi
      # assuming everything is sane, and the input points to a ccd fixed stack, 
      # we can do a few steps automatically.
      ln -s ../../fixedStacks/${inp}.fixed.rot ${pName}.st
      prev_Iter=$((${iter}-1))
      prev_wDir="round_${prev_Iter}"


      # copy the previous tilt angles to rawtlt
      cp ../${prev_wDir}/INP-${prev_Iter}.tlt ./${pName}.rawtlt


      prev_Iter2=$((${iter}-2)) 
      prev_wDir2="round_${prev_Iter2}"

      if [[ $iter -eq 2 ]] ; then
        cp ../${prev_wDir}/INP-${prev_Iter}.inpXF ${pName}.inpXF
      else 
        # Combine previous results
        xfproduct -InputFile1 ../${prev_wDir}/INP-${prev_Iter}.inpXF \
                  -InputFile2 ../${prev_wDir}/INP-${prev_Iter}.tltxf \
                  -OutputFile ${pName}.inpXF
      fi

    fi # if condition on first iteration or not
      

    newstack \
      -InputFile	${pName}.st \
      -OutputFile ${pName}.preali \
      -TransformFile ${pName}.inpXF \
      -BinByFactor ${iBin} \
      -AntialiasFilter 5 \
      -FloatDensities 2
       
    nPrjs=$(wc ${pName}.rawtlt | awk '{print $1}')

    rm -f ${pName}.prexg
    for iPrj in $(seq 1 $nPrjs) ; do
      echo "1.0 0.0 0.0 1.0 0.0 0.0" >> ${pName}.prexg
    done 


    {
     echo '#!/bin/bash'
     echo 'tiltxcorr \'
     echo "-OverlapOfPatchesXandY	${ptOverlap},${ptOverlap} \ "
     echo "-IterateCorrelations	1 \ "
     echo "-ImagesAreBinned	${iBin} \ "
     echo "-InputFile	${pName}.preali \ "
     echo "-OutputFile	${pName}_pt.fid \ "
     echo "-PrealignmentTransformFile	${pName}.prexg \ "
     echo "-TiltFile	${pName}.rawtlt \ "
     echo "-FilterSigma1	${hpCut} \ "
     echo "-FilterRadius2	${lpCut} \ "
     echo "-FilterSigma2	0.05 \ "
     echo "-ShiftLimitsXandY ${limitShifts},${limitShifts} \ "
     echo "-PadsInXandY	128,128 \ "
     echo "-SizeOfPatchesXandY	${ptSize},${ptSize} \ "
     echo "-BordersInXandY ${ptBorder},${ptBorder} \ "
     echo "-CorrelationCoefficient \ "
     echo "-RotationAngle	0.0 \ "
     echo "-AngleOffset ${tiltAngleOffset} "
    } > ./.${inp}_xcorr.sh
    # For some reason specifying BordersInXorY or RotationAngle would terminate
    # input. Get rid of the former, and put rotation at the end.
    tiltxcorr \
      -OverlapOfPatchesXandY	${ptOverlap},${ptOverlap} \
      -IterateCorrelations	1 \
      -ImagesAreBinned	${iBin} \
      -InputFile	${pName}.preali \
      -OutputFile	${pName}_pt.fid \
      -PrealignmentTransformFile	${pName}.prexg \
      -TiltFile	${pName}.rawtlt \
      -FilterSigma1	${hpCut} \
      -FilterRadius2	${lpCut} \
      -FilterSigma2	0.05 \
      -ShiftLimitsXandY ${limitShifts},${limitShifts} \
      -PadsInXandY	128,128 \
      -SizeOfPatchesXandY	${ptSize},${ptSize} \
      -BordersInXandY ${ptBorder},${ptBorder} \
      -CorrelationCoefficient \
      -RotationAngle	0.0 \
      -AngleOffset ${tiltAngleOffset}
      
   
      # This could be specified with "-LengthAndOverlap" in tiltxorr
     CONTOUR_LENGTH="-LengthOfPieces 5"
     imodchopconts \
      -InputModel ${pName}_pt.fid \
      -OutputModel ${pName}.fid \
      -NumberOfPieces 2

    # finally run the alignment
    # Most of these are stock parameters. Note the rotation is set to 0 since
    # the stack is resampled. Increase robust downweighting k --> 0.5
    # Don't fit local tilts, and for the test set which has 100 prjs, increase
    # the default tilt group from 5 --> 8 to address the low angle switch.
    if [[ $iter -gt 2 ]] ; then

      # Run with local
      tiltalign \
        -ModelFile	${pName}.fid \
        -ImageFile	${pName}.preali \
        -ImagesAreBinned	${iBin} \
        -OutputModelFile	${pName}.3dmod \
        -OutputResidualFile	${pName}.resid \
        -OutputFidXYZFile	${pName}fid.xyz \
        -OutputTiltFile	${pName}.tlt \
        -OutputTransformFile	${pName}.tltxf_nonScaled \
        -RotationAngle	0.0 \
        -TiltFile	${pName}.rawtlt \
        -AngleOffset ${tiltAngleOffset} \
        -RotOption	1 \
        -RotDefaultGrouping	3 \
        -TiltOption	${TILT_OPTION} \
        -TiltDefaultGrouping	3 \
        -MagOption	${MAG_OPTION} \
        -MagDefaultGrouping	3 \
        -BeamTiltOption	0 \
        -ResidualReportCriterion	1.0 \
        -SurfacesToAnalyze	0 \
        -RobustFitting	1\
        -MetroFactor	0.25 \
        -MaximumCycles	1000 \
        -KFactorScaling	1.0 \
        -NoSeparateTiltGroups	2 \
        -AxisZShift	1000 \
        -ShiftZFromOriginal      1 \
        -LocalAlignments	\
        -OutputLocalFile	${pName}_local.xf \
        -TargetPatchSizeXandY	${ptSize},${ptSize} \
        -MinSizeOrOverlapXandY	0.5,0.5 \
        -MinFidsTotalAndEachSurface	8,3 \
        -FixXYZCoordinates	1 \
        -LocalRotOption	1 \
        -LocalRotDefaultGrouping	6 \
        -LocalTiltOption	0 \
        -LocalTiltDefaultGrouping	6 \
        -LocalMagOption	1 \
        -LocalMagDefaultGrouping	7  \
        -AngleOffset ${tiltAngleOffset} > tiltAlign.log
    else
      # run without local until an intial stable global solution is found
      tiltalign \
          -ModelFile	${pName}.fid \
          -ImageFile	${pName}.preali \
          -ImagesAreBinned	${iBin} \
          -OutputModelFile	${pName}.3dmod \
          -OutputResidualFile	${pName}.resid \
          -OutputFidXYZFile	${pName}fid.xyz \
          -OutputTiltFile	${pName}.tlt \
          -OutputTransformFile	${pName}.tltxf_nonScaled \
          -RotationAngle	0.00 \
          -TiltFile	${pName}.rawtlt \
          -AngleOffset	${tiltAngleOffset} \
          -RotOption	1 \
          -RotDefaultGrouping	3 \
          -TiltOption	${TILT_OPTION} \
          -TiltDefaultGrouping	3 \
          -MagOption	${MAG_OPTION} \
          -MagDefaultGrouping	3 \
          -BeamTiltOption	0 \
          -ResidualReportCriterion	1.0 \
          -SurfacesToAnalyze	0 \
          -RobustFitting	1\
          -MetroFactor	0.25 \
          -MaximumCycles	1000 \
          -KFactorScaling	1.0 \
          -NoSeparateTiltGroups	2 \
          -AxisZShift	1000 \
          -ShiftZFromOriginal 1  \
          -AngleOffset ${tiltAngleOffset} > tiltAlign.log
     fi

    # This shifts could be scaled using xfproduct, but just keep this simple.
    # Local shifts have scaling info in header, so no worries there.
    awk -v BINVAL=${iBin} '{print $1,$2,$3,$4,BINVAL*$5,BINVAL*$6}' ${pName}.tltxf_nonScaled > ${pName}.tltxf

    iter=$(($iter+1))
    cd ..

  done # bin iterations
done # bin loop
  


finalRound=$(($iter-1))
for j in $(seq 2 $finalRound) ; do
  rm round_${j}/INP-${j}.preali
done


echo $inp >> ../fixedStacks/baseName.list 
##if [[ $FORCE_HIGH_PASS ]] ; then
##  rm ../fixedStacks/${inp}.fixedNP.rot 
##fi


cp round_${finalRound}/INP-${finalRound}_local.xf ../fixedStacks/${inp}.local
xfproduct -InputFile1  round_${finalRound}/INP-${finalRound}.inpXF -InputFile2 round_${finalRound}/INP-${finalRound}.tltxf -OutputFile withOutRot.xf
# Now add in the prerotation
xfproduct -InputFile1 preRotXf.xf -InputFile2 withOutRot.xf -OutputFile ../fixedStacks/${inp}.xf
cp round_${finalRound}/INP-${finalRound}.tlt ../fixedStacks/${inp}.tlt

cd ..




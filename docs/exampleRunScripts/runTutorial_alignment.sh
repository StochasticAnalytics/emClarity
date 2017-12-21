skipThis=1

# If you have fewer physical cores available than requested, your distro will complain.
nCores=24


### If you have run ctf estimation, then you could skip this as the first 
### aligned stacks are alread created.
# emClarity ctf update param0.m emClarity_tutorial${iTomo} full

if [[ $skipThis -eq 0 ]]; then
exit

fi # end of skip, if you want to skip a section, just move this line to that point.

  emClarity init param0.m
    
  emClarity ctf 3d param0.m 
    
  for i in 0 1 2; do
    if [[ $i -eq 0 ]] ; then
      ST='NoAlignment'
    else
      ST='RawAlignment'
    fi

    if [[ $i -ge 0 ]] ; then   emClarity avg param${i}.m ${i} $ST ; fi
      
    if [[ $i -ge 0 ]] ; then  emClarity alignRaw param${i}.m ${i}; fi
      
  done

   
   emClarity tomoCPR param2.m 2 4
  
   # We need to re-create the aligned stacks with the updated alignments from tomoCPR
   emClarity ctf update param3.m -1 full

   emClarity ctf 3d param3.m

  
  # Do the same thing, but now at a binning of 3, adding a duplicate removal step.
  for i in 3 4 5; do
    if [[ $i -ge 3 ]] ; then  emClarity avg param${i}.m ${i} RawAlignment ; fi
    
    if [[ $i -ge 3 ]] ; then  emClarity alignRaw param${i}.m ${i} ; fi
    
  done

  # The ribosomes don't tend to drift too much, but make sure we don't have 
  # any duplicate sub-tomogramnms
  emClarity removeDuplicates param5.m 5 

  emClarity tomoCPR param5.m 5 4
  
  emClarity ctf update param6.m -1 full   
  
  emClarity ctf 3d param6.m

  

  # Again but now at a binning of 2
  for i in 6 7 8 ; do
    if [[ $i -ge 6 ]] ; then  emClarity avg param${i}.m ${i} RawAlignment ; fi
    
    if [[ $i -ge 6 ]] ; then  emClarity alignRaw param${i}.m ${i} ; fi
    
  done
 
  emClarity tomoCPR param8.m 8 4
  
  emClarity ctf update param9.m -1 full  
  
  emClarity ctf 3d param9.m 

  
  # Finally, go to full sampling

  for i in 9 10 11 12  ; do

    if [[ $i -ge 9 ]] ; then  emClarity avg param${i}.m ${i} RawAlignment ; fi
    
    if [[ $i -ge 9 ]] ; then  emClarity alignRaw param${i}.m ${i} ; fi
    
  done

  # A regular average to get a final estimate of the SSNR
  emClarity avg param13.m 13 RawAlignment
  # Align and merge the two half-sets applying requested b-factors.
  emClarity avg param13.m 13 FinalAlignment


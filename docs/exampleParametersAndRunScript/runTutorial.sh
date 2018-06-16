skipThis=1
runThis=1

#
#



if [[ ${skipThis} -eq 0 ]] ; then
  exit

fi #
   emClarity init param0.m ; [[ $? -ne 0 ]] && exit

   emClarity ctf update param0.m ; [[ $? -ne 0 ]] && exit

   emClarity ctf 3d param0.m ; [[ $? -ne 0 ]] && exit
  


  for i in 0 1 2 ; do
    if [[ $i -eq 0 ]] ; then
      ST='NoAlignment'
    else
     ST='RawAlignment'
  fi

    if [[ $i -ge 0 ]] ; then   emClarity avg param${i}.m ${i} $ST ; fi; [[ $? -ne 0 ]] && exit

    if [[ $i -ge 0 ]] ; then  emClarity alignRaw param${i}.m ${i}; fi; [[ $? -ne 0 ]] && exit

  done



  emClarity removeDuplicates param${i}.m ${i} ; [[ $? -ne 0 ]] && exit

   emClarity tomoCPR param${i}.m ${i}  ; [[ $? -ne 0 ]] && exit

   emClarity ctf update param$((${i}+1)).m  ; [[ $? -ne 0 ]] && exit

   emClarity ctf 3d param$((${i}+1)).m; [[ $? -ne 0 ]] && exit
  


  for i in 3 4 5 ; do

    if [[ $i -ge 3 ]] ; then  emClarity avg param${i}.m ${i} RawAlignment ; fi; [[ $? -ne 0 ]] && exit
    
    if [[ $i -ge 3 ]] ; then  emClarity alignRaw param${i}.m ${i} ; fi; [[ $? -ne 0 ]] && exit

  done



   emClarity removeDuplicates param${i}.m ${i} ; [[ $? -ne 0 ]] && exit

   emClarity tomoCPR param${i}.m ${i}  ; [[ $? -ne 0 ]] && exit

   emClarity ctf update param$((${i}+1)).m ; [[ $? -ne 0 ]] && exit

   emClarity ctf 3d param$((${i}+1)).m; [[ $? -ne 0 ]] && exit




  for i in 6 7 8 ; do
    if [[ $i -ge 6 ]] ; then  emClarity avg param${i}.m ${i} RawAlignment ; fi; [[ $? -ne 0 ]] && exit
    
    if [[ $i -ge 6 ]] ; then  emClarity alignRaw param${i}.m ${i} ; fi; [[ $? -ne 0 ]] && exit
    
  done

   emClarity removeDuplicates param${i}.m ${i} ; [[ $? -ne 0 ]] && exit

   emClarity tomoCPR param${i}.m ${i}  ; [[ $? -ne 0 ]] && exit

   emClarity ctf update param$((${i}+1)).m ; [[ $? -ne 0 ]] && exit

   emClarity ctf 3d param$((${i}+1)).m ; [[ $? -ne 0 ]] && exit



  for i in 9 10 11 12 ; do
    
    if [[ $i -ge  9 ]] ; then  emClarity avg param${i}.m ${i} RawAlignment ; fi; [[ $? -ne 0 ]] && exit

    if [[ $i -ge  9 ]] ; then  emClarity alignRaw param${i}.m ${i} ; fi; [[ $? -ne 0 ]] && exit

  done


  emClarity avg param13.m 13 RawAlignment; [[ $? -ne 0 ]] && exit
  emClarity avg param13.m 13 FinalAlignment; [[ $? -ne 0 ]] && exit
 # fi # end of skipThis
#fi # end of runTHis

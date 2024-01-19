function [ pool ] = EMC_parpool(nWorkers)
%Make a local copy of the default cluster, and modify the job storage
%location


local_cache_root = getenv('MCR_CACHE_ROOT');
if isempty(local_cache_root)
  fprintf('\nWARNING: the variable MCR_CACHE_ROOT is not set!\n');
  EMC_ROOT = getenv('EMC_CACHE_DIR');
  if isempty(EMC_ROOT)
    error('\nWARNING: the variable EMC_CACHE_ROOT is not set!\n');
  else
    error('\The variable EMC_CACHE_ROOT %s\n',EMC_ROOT);
  end
end


[~,emc_rand_name,~] = fileparts( local_cache_root );


% check to see if we have already made a parpool in this call to emClarity.
profile_does_not_exist = true;
current_profiles = parallel.clusterProfiles;
for iProf = 1:length(current_profiles)
  if strcmp(current_profiles{iProf},emc_rand_name)
    profile_does_not_exist = false;
  end
end

% if the profile doesn't exist, create it
if (profile_does_not_exist)
  emc_parcluster = parcluster(parallel.defaultClusterProfile);
  emc_parcluster.JobStorageLocation = local_cache_root;
  saveAsProfile(emc_parcluster, emc_rand_name);
end

[ pool ] = parpool(emc_rand_name, nWorkers);

end


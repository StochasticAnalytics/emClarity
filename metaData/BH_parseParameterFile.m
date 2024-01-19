function [ emc ] = BH_parseParameterFile( PARAMETER_FILE )
%Parse a parameter file & check for valid parameters.
%   experimental

fileID = fopen(PARAMETER_FILE,'r');

p = textscan(fileID, '%s', 'CommentStyle',{'%'},'Delimiter','\n', ...
  'TreatAsEmpty',{' '});

nParam = 1;
p2 = cell(1,1);
% Get rid of empty lines
for i = 1:size(p{1},1)
  if ~isempty(p{1}{i})
    p2{nParam,1} = p{1}{i};
    nParam = nParam + 1;
  end
end
clear p

emc = struct();
% Check that all paramters are name: value pairs
stringValues = {'subTomoMeta'; ...
  'Ali_mType';'Cls_mType';'Cls_mType';'Raw_mType';'Fsc_mType'; ...
  'Pca_distMeasure';'Kms_mType';'flgPrecision';'Tmp_xcfScale';...
  'fastScratchDisk';'Tmp_eraseMaskType';'startingDirection';'Peak_mType';'symmetry'};
for i = 1:size(p2,1)
  pNameVal = strsplit(p2{i,1},'=');
  if length(pNameVal) == 1
    error('Could not split Name=Value pair for\n\t %s',char(pNameVal))
  elseif length(pNameVal) > 2
    error('To many colons in\n\t %s',char(pNameVal))
  else
    
    if any(strcmp(stringValues, pNameVal{1}))
      emc.(pNameVal{1}) = pNameVal{2};
    else
      emc.(pNameVal{1}) = EMC_str2double(pNameVal{2});
    end
  end
end

% Now check for optional parameters

% Early development parameter, used to store more than one orientation during template matching
% and use for further refinement.
if isfield(emc, 'nPeaks')
  EMC_assert_integer(emc.nPeaks, 1);
else
  emc.('nPeaks') = 1;
end

% Used when cutting out subtomos for further processing. Adds extra padding to anticipate shifts etc.
% This has not been well tested

% When used in average3d, this value is stored in the subTomoMeta.
if isfield(emc, 'CUTPADDING')
  EMC_assert_integer(emc.CUTPADDING, 1);
else
  emc.('CUTPADDING') = 20;
end

if isfield(emc, 'whitenPS')
  EMC_assert_numeric(emc.whitenPS, 3)
  emc.('wiener_constant') = emc.whitenPS(3);
else
  emc.('whitenPS') = [0.0,0.0,0.0];
  emc.('wiener_constant') = 0.0;
end

% Default bfactor applied to the re-weighting when generating the fully corrected volumes.
% positive corresponds to a sharpening, negative to a low-pass.
if isfield(emc, 'Fsc_bfactor')
  EMC_assert_numeric(emc.Fsc_bfactor)
else
  emc.('Fsc_bfactor') = 40.0;
end

% Used to downweight higher frequencies based on relative CCC scores.
% Based on one of Niko's papers, but catching some edge cases for tomo.
% Overwritten if cycle == 0 as the scores from template matching do not work for this metric as they are SNR not CCC.
% TODO: get rid of the flg prefix
if isfield(emc, 'flgQualityWeight')
  EMC_assert_numeric(emc.flgQualityWeight, 1)
else
  emc.('flgQualityWeight') = 5.0;
end

% Experimental downweighting of higher frequency info farther from focus.
% Could also consider filtering pre reconstruction
% Filtering by defocus using exp[-(%d*(argmax(def-1,0,5).*q)^%d)]\n',flgFilterDefocus);
if isfield(emc,'filterDefocus')
  EMC_assert_numeric(emc.filterDefocus, 2)
else
  emc.filterDefocus = [0.0, 0.0];
end

if isfield(emc,'flgCutOutVolumes')
  EMC_assert_boolean(emc.flgCutOutVolumes)
else
  emc.flgCutOutVolumes = false;
end


if isfield(emc,'track_stats')
  EMC_assert_boolean(emc.track_stats)
else
  emc.track_stats = false;
end


% Check and override the rotational convention to get helical averaging.
% Replaces the former hack of adding a fifth dummy value to the angular search
if isfield(emc,'doHelical')
  EMC_assert_boolean(emc.doHelical)
else
 emc.doHelical = false;
end

end


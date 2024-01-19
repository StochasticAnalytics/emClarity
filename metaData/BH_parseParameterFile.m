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

if ~isfield(emc, 'nPeaks')
  emc.('nPeaks') = 1;
end

end


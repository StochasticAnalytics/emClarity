function [ pStruct ] = BH_parseParameterFile( PARAMETER_FILE )
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

pStruct = struct();
% Check that all paramters are name: value pairs
stringValues = {'subTomoMeta'; ...
                'Ali_mType';'Cls_mType';'Cls_mType';'Raw_mType';'Fsc_mType'; ...
                'Pca_distMeasure';'Kms_mType';'flgPrecision';'Tmp_xcfScale';...
                'fastScratchDisk';'Tmp_eraseMaskType'};
for i = 1:size(p2,1)
  pNameVal = strsplit(p2{i,1},'=');
  if length(pNameVal) == 1
      error('Could not split Name=Value pair for\n\t %s',char(pNameVal))
  elseif length(pNameVal) > 2
    error('To many colons in\n\t %s',char(pNameVal))
  else

    if any(strcmp(stringValues, pNameVal{1}))
      pStruct.(pNameVal{1}) = pNameVal{2};
    else
      pStruct.(pNameVal{1}) = str2num(pNameVal{2});
    end
  end
end


end


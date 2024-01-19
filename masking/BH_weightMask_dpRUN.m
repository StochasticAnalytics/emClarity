function [] = BH_weightMask_dpRUN(PARAMETER_FILE,CYCLE,PRFX_OUT,SYMMETRY,varargin)
%UNTITLED Summary of this function goes heres
%   Detailed explanation goes here

ctfScaleFactor = 1;
tomoStart=0;
tomoStop = 0;
CYCLE = EMC_str2double(CYCLE);
SYMMETRY = EMC_str2double(SYMMETRY);
nargin
if nargin > 4
  varargin{1} = EMC_str2double(varargin{1})
  if length(varargin{1}) > 1
    tiltStart = varargin{1}(2); %ctfScaleFactor = varargin{1}(2)\
    tiltStop = varargin{1}(3);
  else 
    tiltStart = 0;
    tiltStop = 0;
  end
  useGPU = varargin{1}(1);
else 
  useGPU = -1;
end
tiltStart
tiltStop
gpuIDX = BH_multi_checkGPU(useGPU);
gDev = gpuDevice(gpuIDX);

emc = BH_parseParameterFile(PARAMETER_FILE);
cycleNumber = sprintf('cycle%0.3d',CYCLE);

load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');


geom = subTomoMeta.(cycleNumber).Avg_geometry;
% Keep this in to include optional run of applying the weights
fVal = subTomoMeta.(cycleNumber).fitFSC.Raw1;
mVal = subTomoMeta.(cycleNumber).fitFSC.MaskRaw1;
aVal = subTomoMeta.(cycleNumber).fitFSC.ResampleRaw1;

wgtSize = subTomoMeta.(cycleNumber).class_0_Locations_REF_ODD_Wgt{2}{1}(2:2:6);
samplingRate = emc.('Ali_samplingRate');

tomoList = fieldnames(geom);
if ( tiltStart )
  tomoList = tomoList(tiltStart:tiltStop);
  fprintf('truncate tomolist to %d 5d\n',tomoStart,tomoStop);
end

[ctfWeights] = BH_weightMask_dp(subTomoMeta,wgtSize,samplingRate, ...
                                {tomoList,geom},'double','GPU', ...
                                ctfScaleFactor);
                              
if SYMMETRY > 1
  for iGold = 1:2
    ctfWeights{iGold} = gather(BH_resample3d(ctfWeights{iGold},[0,0,0], ...
                                            [0,0,0],{'Bah',SYMMETRY,'linear',1}, ...
                                            'GPU','forward'));
  end
end

for iGold = 1:2
  SAVE_IMG(MRCImage(ctfWeights{iGold}),sprintf('%s_%d_%d_wgt_%d.mrc',PRFX_OUT,tiltStart,tiltStop,iGold));
end

clear ctfWeights
end


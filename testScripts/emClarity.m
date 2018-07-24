function [ varargout ] = emClarity( varargin )
%Wrapper to run BH_subtomo programs.
%   Detailed explanation goes here

% Disable warnings
warning off
cudaStart='';
cmdIN = sprintf('emClarity %s ',varargin{1});

% first argument is the program to run or "help" to print a list of available
% options.
setenv('MATLAB_SHELL','/bin/bash');
[sysERR] = system('mkdir -p logFile');
if ( sysERR )
  fprintf('system command for mkdir failed for some reason\n');
  unix('mkdir -p logFile');
end
diary('logFile/emClarity.logfile')
timeStart = datetime();
fprintf('\n\t\t***************************************\n\n');
fprintf('emClarity version %s\n', varargin{1});
fprintf('run starting on %s\n', timeStart);
fprintf('cmd ');
for iArg = 2:length(varargin)
  cmdIN = [cmdIN sprintf('%s ',varargin{iArg})];
  fprintf('%s ', varargin{iArg});
end
fprintf('\n\n')
fprintf('\t\t***************************************\n\n');
% Get rid of the shorthead passed in by the emClarity script.
varargin = varargin(2:end);
if strcmpi(varargin{1},'h')
  varargin{1} = 'help';
end
% nargin doesn't adjust for the shift.
nArgs = length(varargin);
% For some reason, this doesn't fall through, so error if < arg
if nArgs > 1
  if ~strcmp(varargin{2},'check')
    checkHelp = ~strcmpi(varargin{2},'help') || ~strcmp(varargin{2},'experimental');
  end
elseif nArgs == 0
  error('\n\n\tRun with help for a list of functions\n\n');
%   checkHelp = 0;
end

multiGPUs = 1;
if nArgs > 1 && checkHelp
  switch varargin{1}
    case 'ctf'
      pBH = emC_testParse(varargin{3});
    case 'rescale'
      % nothing to parse
      multiGPUs = 0;
    case 'alignFrames'
      error('alignFrames is not yet in production. Soon though!')
      % nothing to parse
      multiGPUs =0;
    case 'benchmark'
      % nothing to parse
      multiGPUs= 0;
    case 'removeNeighbors'
      % nothing to parse
      multiGPUs = 0;

    otherwise
      pBH = emC_testParse(varargin{2});
  end
  if ( multiGPUs )
    % wanted num gpus
    nGPUs_wanted = pBH.('nGPUs');
    cudaStart = getenv('CUDA_VISIBLE_DEVICES')
    nGPUs_visible = gpuDeviceCount;
    if nGPUs_visible < nGPUs_wanted
      error('\n\n\t\t%d gpus requested but only %d are visible to the system\n\n',nGPUs_wanted, nGPUs_visible);
    elseif nGPUs_visible > nGPUs_wanted
      % Select largest mem visible
      fprintf('\nThere are more gpus visible than requested, selecting the largest memory devices\n');
      select_gpus(nGPUs_wanted,nGPUs_visible,cmdIN);
    else
     fprintf('\nThe number of gpus requested matches those visible to the system\n');
    end
  end
end

switch varargin{1}
  case 'help'
    fprintf(['\nAvailable commands (case sensitive):\n\n',...
             '\ncheck - system check for dependencies\n',...
             '\ninit - create a new project from template matching results.\n',...
             '\navg - average subtomograms\n',...
             '\nfsc - calculate the fsc\n',...
             '\nalignRaw - align one or more references against individual subtomograms.\n',...
             '\npca - reduce dimensionality prior to clustering, possibly on smaller subset of data.\n',...
             '\ncluster - use one of a number of approaches to sort populations.\n',...
             '\nskip - after avering classes & possible removing some, skip to next cycle.\n',...
             '\ngeometry - edit or analyze the experimental metadata.\n',...
             '\nctf - estimate, correct, or refine the CTF.\n',...
             '\ntomoCPR - tomogram constrained projection refinement\n',...
             '\ntemplateSearch - template matching/ global search\n',...
             '\ncleanTemplateSearch - clean search results based on neighbor constraints\n',...
             '\nrescale - change the mag on a volume\n',...
             '\nremoveDuplicates - remove subtomos that have migrated to the same position\n',...
             '\nremoveNeighbors - clean templateSearch results based on lattice constraints\n']);
           
% Currently disabled options. Multi-reference alignment 
% % %                       '\nalignRef - align one or more references against a primary ref\n',...
% % %              '             optionally add multiple instances aligned to each\n',...
% % %              '             to its respective reference.\n',...
% % %              '\nalignCls - align refs from alignRef to a usually much larger number of class averages.\n',...
  case 'experimental'
    print_experimental_options();
  case 'check'
    if nArgs > 1
      fprintf('check takes no arguments');
    else
      BH_checkInstall(  )
    end
  case 'init'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 2 && length(varargin) ~= 3
      fprintf(['\nparam.m [tomoCpr iter, for continuing after second globalsearch]\n']);
    elseif length(varargin) == 3
      emC_testParse(varargin{2})
      BH_geometryInitialize(varargin{2},varargin{3});
    else
      emC_testParse(varargin{2})
      BH_geometryInitialize(varargin{2});
    end   
  case 'removeNeighbors'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 6 
      fprintf(['\npixelSize CYCLE distanceCutoff (Ang) angleCutoff (Deg) N-neighbors\n']);
    else
      %emC_testParse(varargin{2})
      BH_geometry_Constraints(varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
    end    
    
    
  case 'skip'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 3
      fprintf(['\nparam.m iter\n']);
    else
      emC_testParse(varargin{2})    
      BH_skipClassAlignment(varargin{2},varargin{3},'ClassAlignment','1');
      BH_skipClassAlignment(varargin{2},varargin{3},'RawAlignment','1');
    end
  case 'rescale'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 6
      fprintf(['\nfileNameIN fileNameOut angPixIN angPixOut cpu/GPU\n']);
    else
      mag = str2double(varargin{4})/str2double(varargin{5});
      BH_reScale3d(varargin{2},varargin{3},mag,varargin{6});
    end
  case 'benchmark'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 4
      fprintf(['\nfileNameOut fastScratchDisk nWorkers']);
    else
      BH_benchmark(varargin{2},varargin{3},varargin{4});
    end    
  case 'calcWeights'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 6
      fprintf('%f\n',length(varargin))
      fprintf(['\nparam.m cycle prefixOUT symmetry [gpuIDX, tiltStart, tiltStop]\n']);
    else
      
      BH_weightMask_dpRUN(varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
    end
  case 'avg'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 4
      fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'stage of alignment\n',...
               '  raw (post raw alignment)\n',...
               '  cluster_ref or cluster_cls (post classification)\n']);
    else
        emC_testParse(varargin{2})
        BH_average3d(varargin{2}, varargin{3}, varargin{4});
    end 
  case 'fsc'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       (length(varargin) ~= 4 &&  length(varargin) ~= 6)
      fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'stage of alignment\n',...
               '  raw (post raw alignment)\n',...
               '  cluster_ref or cluster_cls (post classification)\n']);  
    elseif length(varargin) == 4
      emC_testParse(varargin{2})
      BH_fscGold_class(varargin{2}, varargin{3}, varargin{4});
    else
      BH_fscGold_class(varargin{2}, varargin{3}, varargin{4},varargin{5},varargin{6});
    end
  case 'alignRaw'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       (length(varargin) ~= 3 && length(varargin) ~= 4)
     fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'stage of alignment\n']);
    else
        emC_testParse(varargin{2})
        if length(varargin) == 3
          BH_alignRaw3d(varargin{2}, varargin{3});
          
        else
          BH_alignRaw3d(varargin{2}, varargin{3},varargin{4});
        end
    end
    
  case 'alignRef'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 3
     fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'stage of alignment\n']);
    else
        emC_testParse(varargin{2})
        BH_alignReferences3d(varargin{2}, varargin{3});
    end    
  case 'alignCls'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 3
     fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'stage of alignment\n']);
    else
        emC_testParse(varargin{2})
        BH_alignClassRotAvg3d(varargin{2}, varargin{3});
    end    
%   case 'alignFrames'
%     if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
%        (length(varargin) < 9 && length(varargin) > 11)
%      fprintf(['\nnameIN\n',...
%                'nameOUT\n',...
%                'gpuIDX\n',...
%                'FPN\n',...
%                'pixelSizeIN\n',...
%                'pixelSizeOUT\n',...
%                'overSampleBy\n',...
%                'doLocal [particleRadiusAng, maxRes]\n']);
%     else
%       switch length(varargin)
%         case 9
%           BH_alignSubFramesTot(varargin{2}, varargin{3}, varargin{4},...
%                               str2num(varargin{5}),...
%                               str2num(varargin{6}),...
%                               str2num(varargin{7}),...
%                               str2num(varargin{8}),...
%                               str2num(varargin{9}));
%         case 10
%           BH_alignSubFramesTot(varargin{2}, varargin{3}, varargin{4},...
%                               str2num(varargin{5}),...
%                               str2num(varargin{6}),...
%                               str2num(varargin{7}),...
%                               str2num(varargin{8}),...
%                               str2num(varargin{9}),...
%                               str2num(varargin{10}));    
%         case 11
%           BH_alignSubFramesTot(varargin{2}, varargin{3}, varargin{4},...
%                               str2num(varargin{5}),...
%                               str2num(varargin{6}),...
%                               str2num(varargin{7}),...
%                               str2num(varargin{8}),...
%                               str2num(varargin{9}),...
%                               str2num(varargin{10}),...
%                               str2num(varargin{11}));
%       end 
%     end      
  case 'pca'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 4
     fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'randomSubset\n']);
             %  'use focused mask\n',...
             %  '  1 from standard devation\n',...
             %  '  2 from variance\n',...
             %  '  3 user supplied (not recommended)\n']);
    else
        emC_testParse(varargin{2})

        maskVal = 0; %str2double(varargin{5});
        if (maskVal == 3)
          % run on full or random subset using user supplied mask
          BH_pcaPub(varargin{2}, varargin{3}, '-3')
        else
          % run on full or randomsubset with no variance mask availble
          BH_pcaPub(varargin{2}, varargin{3}, '0')
          if (maskVal)
            % re-run on full or randomsubset now using variance or stddev mask
            BH_pcaPub(varargin{2}, varargin{3}, sprintf('%d',-1.*maskVal))
          end
        end

        if (str2double(varargin{4}))
          % project onto full set
          BH_pcaPub(varargin{2}, varargin{3}, '1')
        end       
       
    end    
  case 'cluster'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 3
     fprintf(['\nparam.m\n',...
               'cycle number\n']);
    else
        emC_testParse(varargin{2})
        BH_clusterPub(varargin{2}, varargin{3});
    end      
  case 'ctf'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') 
      fprintf(['\nestimate\n',...
               '  param.m tiltBaseName\n',...
               '\nrefine\n',...
               '  param.m tiltBaseName gpuIDX\n',...
               '\nupdate\n',...
               '  param.m tiltBaseName (full,refine,update)\n',...
               '\ncorrect\n',...
               '  param.m precision      usuable-area nWorkers\n',...
               '         (single,double) [nx,ny,nz] \n',...
               '\n3d\n',...
               '  param.m [/local/Scratch]\n']);
    else
      
      emC_testParse(varargin{3})    

      switch varargin{2}
        case 'estimate'
          if nargin == 5
            BH_ctf_Estimate(varargin{3},varargin{4});
          else
            BH_ctf_Estimate(varargin{3},varargin{4},varargin{5});
          end
        case 'refine'
          if length(varargin) == 4
            error('You need to specify a gpuIDX')
          end          
          BH_ctf_Refine2(varargin{3},varargin{4},varargin{5});
        case 'update'
          if length(varargin) > 3
            error('\n\nYou now only need to specify %s parameter file.\n\n','the')
          end          
          BH_ctf_Updatefft(varargin{3},'-1','full');
        case 'correct'
          BH_ctf_Correct(varargin{3},varargin{4},varargin{5},varargin{6},varargin{7});
        case '3d'
          if nArgs == 5
            % Not a public option, start from tilt # (of nTilts)       
            BH_ctf_Correct3d(varargin{3},varargin{4},varargin{5});            
          elseif nArgs == 4
            BH_ctf_Correct3d(varargin{3},varargin{4});
          else
            BH_ctf_Correct3d(varargin{3});
          end
            
        otherwise
          error('ctf operations are estimate,refine,update,correct, or 3d.');
      end
    end    
  case 'tomoCPR'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       ( length(varargin) < 3 || length(varargin) > 4 )
      fprintf(['\nparam.m\n',...
               'cycle number\n',...                                  
               'nTiltStart\n']);
    else
        emC_testParse(varargin{2})
        if length(varargin) == 4
          tiltStart = str2double(varargin{4});
        else
          tiltStart = 1;
        end
        BH_synthetic_mapBack(varargin{2}, varargin{3},tiltStart);
    end        
  case 'removeDuplicates'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 3
      fprintf(['\nparam.m\n',...
               'cycle number\n',...
                ]);
    else
        emC_testParse(varargin{2})
        BH_removeDuplicates(varargin{2}, varargin{3} );
    end       
  case 'geometry'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 7
      fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'stage of alignment\n',...
               'operation []\n',...
               '  SwitchCurrentCycle, UpdateTilts, WriteCsv, RemoveClasses,\n'...
               '  ShiftAll, ShiftBin, ListTomos, RemoveTomos,\n',...
               '  ListPercentiles, RemoveFraction, RandomizeEulers\n',...
               'vectOP [0,0,0]\n',...
               'STD, EVE, ODD\n']);
    else
        emC_testParse(varargin{2})
        BH_geometryAnalysis(varargin{2}, varargin{3},varargin{4}, ...
                            varargin{5}, varargin{6},varargin{7});
    end 
  case 'calcWeights'
    
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       (length(varargin) < 5 && length(varargin) > 7 )
      fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'name prefix\n',...
               'symmetry (in plane) [gpuIDX]']);
    else
      if length(varargin) == 5
         emC_testParse(varargin{2})
         BH_weightMask_dpRUN(varargin{2},str2double(varargin{3}),varargin{4}, ...
                                                     str2double(varargin{5}));
      elseif length(varargin) == 6
         emC_testParse(varargin{2})
         BH_weightMask_dpRUN(varargin{2},str2double(varargin{3}),varargin{4}, ...
                                 str2double(varargin{5}),str2double(varargin{6}));  
      elseif length(varargin) == 7
         emC_testParse(varargin{2})
         BH_weightMask_dpRUN(varargin{2},str2double(varargin{3}),varargin{4}, ...
                                 str2double(varargin{5}), ...
                                 [str2double(varargin{6}),str2double(varargin{7})]);                                 
      end  
    end
    
    
  case 'templateSearch'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       (length(varargin) ~= 7 && length(varargin) ~= 8)
       fprintf(['\nparam.m\n',...
           'tomoName\n',...
           'tomoNumber\n', ...
           'template name\n',...
           'symmetry\n',...
           '[threshold override]\n',...
           'gpuIDX.\n']);
    else
      wedgeType = 4;
      if length(varargin) == 7
        BH_templateSearch3d( varargin{2}, varargin{3},varargin{4}, ...
                           varargin{5}, varargin{6},wedgeType,varargin{7});
      else
        BH_templateSearch3d( varargin{2}, varargin{3},varargin{4}, ...
                           varargin{5}, varargin{6},wedgeType, ...
                           varargin{7},varargin{8});
      end
    end
    
  case 'cleanTemplateSearch'
     if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       length(varargin) ~= 5 
       fprintf(['\npixelSize (Ang)\n',...
           'distance to neightbor (Ang)\n',...
           'angular deviation to neighbor (degrees)\n', ...
           'min number neighbors (one less than expected is usually good)\n']);  
     else
    
      BH_geometry_Constraints(varargin{2}, '0', varargin{3}, varargin{4}, varargin{5});

     end
  otherwise
    error('command --%s-- not recognized. Try "help" for a list.', varargin{1})
end

timeFinish = datetime();
fprintf('\n\t\t***************************************\n\n');
fprintf('run ending on %s\n', timeFinish);
fprintf('\n\n')
fprintf('\t\t***************************************\n\n');


diary off
if (cudaStart)
  setenv('CUDA_VISIBLE_DEVICES',cudaStart);
end
end




function [ pBH  ] = emC_testParse( paramTest )
% Try to parse the parameter file make sure it's okay.
%
% Add some actual error handling here to help trouble shoot.
try
  pBH = BH_parseParameterFile( paramTest );
catch
  error('error parsing parameter file %s\n', paramTest)
end



end

function [] = select_gpus(nGPUs_wanted,nGPUs_visible,cmdIN);

% I don't like this string parsing with system tools
[~,uuidList] = system('nvidia-smi --list-gpus | awk -F "UUID: " ''{print $2}'' | awk -F ")" ''{print $1}''');
uuidlist = strsplit(uuidList);
nGPUs_total = length(uuidlist)-1;
memList = zeros(nGPUs_total,2);
memList(:,1) = 1:nGPUs_total;
for iGPU = 1:nGPUs_total
  memCMD = sprintf('nvidia-smi --id=%s --query-gpu=memory.total --format=csv,noheader,nounits',uuidlist{iGPU});
  [~,memString] = system(memCMD);
  memList(iGPU,2) = str2double(memString);
  fprintf('found gpu with %f Mib memory\n',memList(iGPU,2));
end

memList = sortrows(memList,-2);
devList = '';
for iGPU = 1:nGPUs_total
  if iGPU <= nGPUs_wanted
  devList = strcat(devList,uuidlist{memList(iGPU,1)});
    if iGPU < nGPUs_wanted
      devList = strcat(devList,',');
    elseif iGPU == nGPUs_wanted && nGPUs_total > nGPUs_wanted
      % start list of uuids to explicitly block
      devList = strcat(devList,',-1,');
    end
  elseif iGPU > nGPUs_wanted && nGPUs_total >=nGPUs_wanted
    devList =  strcat(devList,uuidlist{memList(iGPU,1)});
    if iGPU < nGPUs_total
      devList = strcat(devList,',');
    end
  end
end

setenv('CUDA_VISIBLE_DEVICES',char(devList));

fprintf('\n\n\t\tPlease add the following line to your emClarity run script\n\nexport CUDA_VISIBLE_DEVICES=%s\n\n',devList);

exit

end % end select gpus

function print_experimental_options

  fprintf('\n\n\tExperimental Options: use at your own RISK\n');
  fprintf('(\t\tOr better yet, check with ben!\t\t)\n');
  fprintf('\nIf you do use/change any of these, please mention in your methods and EMDB entry!\n');
  fprintf('\n\n----------------------------------\n\n');
  fprintf('\nscaleCalcSize\toversampling of vol for xcorr. Def:\t1.5\n');
  fprintf('\npaddedSize\tpadded size of tiles in ctf estimateion\n');
  
  fprintf('\nflgFscShapeMask\t default 1\n');
  fprintf('\nflgPcaShapeMask\t default 1\n');
  
  fprintf('\nflgQualityWeight\t Downweight high-freq of low scoring sub-tomos. Def:\t4\n');
  fprintf('\ninterpOrder\t Linear interpolation (1) Spline/Fourier (4 - not working, do not use)\n');
  fprintf('\nflgLimitToOneProcess\t For OOM issues in averaging. Boolean Def:\t0\n');
  fprintf('\nflgCenterRefCOM\tShift reference to center of mass. Boolean Def:\t1\n');
  fprintf('\nconserveDiskSpace\n');
  fprintf('\nPca_distMeasure\tMeasure for difference. euclidean, cityblock, correlation, cosine Def:\t sqeuclidean\n');
  fprintf('\nPca_nReplicates\tThe number of times Kmeans is intialized. Def:\t 128\n');
  fprintf('\nflgSymmetrizeSubTomos\tApply symmetry to subtomos in alignment.\nCurrently not worse, not better, but much slower.\n');
  
  fprintf('\ndeltaZTolerance\tallowed defocus variance in ctf estimation:Def:\t100e-9\n');
  fprintf('\nzShift\tselect tiles with a defocus offset. Determine tilt gradient.\n\n');
 

end % end of print experimental options


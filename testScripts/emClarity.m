function [ varargout ] = emClarity( varargin )
%Wrapper to run BH_subtomo programs.
%   Detailed explanation goes here

% Disable warnings
warning off
cudaStart='';
useV2 = false;
cmdIN = sprintf('emClarity %s ',varargin{1});


% first argument is the program to run or "help" to print a list of available
% options.
setenv('MATLAB_SHELL','/bin/bash');
[sysERR] = system('mkdir -p logFile');

if strcmp(varargin{2},'gui')
  emClarityApp;
  return
end
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
varargin;

% For some reason, this doesn't fall through, so error if < arg
notCheckHelp = 1;
if nArgs > 1
  if ~strcmp(varargin{1},'check')
    if strcmp(varargin{1},'v2')
      useV2 = true;
      varargin = varargin(2:end);
    else
      notCheckHelp = ~strcmpi(varargin{2},'help') || ~strcmp(varargin{2},'experimental');
    end
  end
elseif nArgs == 0
  error('\n\n\tRun with help for a list of functions\n\n');
%   checkHelp = 0;
end

multiGPUs = 1;
if nArgs > 1 && notCheckHelp
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
    case 'combineProjects'
      % nothing to parse
      multiGPUs = 0;
    case 'autoAlign'
         
      autoAliPath='/groups/grigorieff/home/himesb/work/emClarity/mexFiles/compiled/emC_autoAlign.sh';
      if isdeployed
        autoAliPath = sprintf('%s%s',ctfroot,autoAliPath);
      end
  
      if ~exist(varargin{3}, 'file')
        fprintf('Did not find your .rawtlt file %s\n',varargin{3});
        error('Expecting tiltName.st tiltName.rawtlt pixelSize (Ang) imageRotation (degrees)');
      end
      if ~exist(varargin{2}, 'file')
        fprintf('Did not find your .st file %s\n',varargin{2});
        error('Expecting tiltName.st tiltName.rawtlt pixelSize (Ang) imageRotation (degrees)');
      end
      
      BH_runAutoAlign(autoAliPath,varargin{2},varargin{3},varargin{4},varargin{5});
           
      return

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
             '\ncombineProjects - combine two or more projects together', ...
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
       length(varargin) < 2 && length(varargin)> 5
      fprintf(['\nparam.m [tomoCpr iter, for continuing after second globalsearch]\n']);
    elseif length(varargin) == 5
      emC_testParse(varargin{2})
      BH_geometryInitialize(varargin{2},varargin{3},varargin{4},varargin{5});
    elseif length(varargin) == 4
      emC_testParse(varargin{2})
      BH_geometryInitialize(varargin{2},varargin{3},varargin{4});
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
  case 'mask'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       (~ismember(length(varargin),[5,8,9]))
     length(varargin)
      fprintf(['\nFor geometric mask:\n', ...
               'fileNameOUT.mrc, pixelSize (Ang), Shape (sphere,cylinder,rectangle), Size/radius/center in pixels: [nX,nY,nZ], [rX,rY,rZ], [cX,cY,cZ], optional: "2d"',...
               '\n\nFor a shape based mask\n', ...
               'fileNameIN.mrc,fileNameOUT.mrc, pixelSize (Ang)\n']);
       
    else
      switch length(varargin) 
        case 5
          maskVol = getVolume(MRCImage(varargin{3}));
          pixelSize = str2double(varargin{5});
          maskVol = BH_mask3d(maskVol,str2double(varargin{5}),'','');
          SAVE_IMG(MRCImage(gather(maskVol)),varargin{4},pixelSize);
        case 8
         pixelSize = str2double(varargin{4});
          maskVol = BH_mask3d(varargin{5},str2num(varargin{6}), ...
                                          str2num(varargin{7}), ...
                                          str2num(varargin{8}));
          SAVE_IMG(MRCImage(gather(maskVol)),varargin{3},pixelSize);
        case 9
           pixelSize = str2double(varargin{4});
          maskVol = BH_mask3d(varargin{5},str2num(varargin{6}), ...
                                          str2num(varargin{7}), ...
                                          str2num(varargin{8}), ...
                                          str2num(varargin{9}));
          SAVE_IMG(MRCImage(gather(maskVol)),varargin{3},pixelSize);
      end
      
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
       length(varargin) ~= 4  && length(varargin) ~= 5
     fprintf(['\nparam.m\n',...
               'cycle number\n',...
               'randomSubset\n']);
             %  'use focused mask\n',...
             %  '  1 from standard devation\n',...
             %  '  2 from variance\n',...
             %  '  3 user supplied (not recommended)\n']);
    else
        emC_testParse(varargin{2})

        if (str2double(varargin{4}))
          % project onto full set
          BH_pcaPub(varargin{2}, varargin{3}, '1');
        else
          
          if length(varargin) == 5
            maskVal = str2double(varargin{5});
          else
            maskVal = 0;
          end

          if (maskVal)
            % re-run on full or randomsubset now using variance or stddev mask
            BH_pcaPub(varargin{2}, varargin{3}, sprintf('%d',-1.*maskVal))
          else
            BH_pcaPub(varargin{2}, varargin{3}, '0')
          end
 
        

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
          if (useV2)
            if nArgs == 4
              BH_ctf_Estimate_2(varargin{3},varargin{4});
            else
              BH_ctf_Estimate_2(varargin{3},varargin{4},varargin{5});
            end            
          else
            if nArgs == 4
              BH_ctf_Estimate(varargin{3},varargin{4});
            else
              BH_ctf_Estimate(varargin{3},varargin{4},varargin{5});
            end
          end
        case 'refine'
          if length(varargin) ~= 4
            error('You need to specify a parameter file and tilt name')
          end          
          BH_ctf_Refine2(varargin{3},varargin{4});
        case 'update'
          if length(varargin) > 3
            error('\n\nYou now only need to specify %s parameter file.\n\n','the')
          end          
          BH_ctf_Updatefft(varargin{3},'-1','full');
        case 'correct'
          BH_ctf_Correct(varargin{3},varargin{4},varargin{5},varargin{6},varargin{7});
        case '3d'
          if nArgs == 6
            % last is a dummy, used for tomoCPR background
            BH_ctf_Correct3d(varargin{3},varargin{4},varargin{5},varargin{6});     
          elseif nArgs == 5
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
  case 'combineProjects'
    BH_combineProjects(varargin{1},varargin(2:end));
    
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
    
    
% % %   case 'templateSearch'
% % %     if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
% % %        (length(varargin) ~= 7 && length(varargin) ~= 8)
% % %        fprintf(['\nparam.m\n',...
% % %            'tomoName\n',...
% % %            'tomoNumber\n', ...
% % %            'template name\n',...
% % %            'symmetry\n',...
% % %            '[threshold override]\n',...
% % %            'gpuIDX.\n']);
% % %     else
% % %       wedgeType = 2;
% % %       if length(varargin) == 7
% % %         BH_templateSearch3d( varargin{2}, varargin{3},varargin{4}, ...
% % %                            varargin{5}, varargin{6},wedgeType,varargin{7});
% % %       else
% % %         BH_templateSearch3d( varargin{2}, varargin{3},varargin{4}, ...
% % %                            varargin{5}, varargin{6},wedgeType, ...
% % %                            varargin{7},varargin{8});
% % %       end
% % %     end
    
  case 'templateSearch'
    if strcmpi(varargin{2},'help') || strcmpi(varargin{2},'h') || ...
       ~ismember(length(varargin),[7,8])
       fprintf(['\nparam.m\n',...
           'tomoName\n',...
           'tomoNumber\n', ...
           'template name\n',...
           'symmetry\n',...
           '[threshold override]\n',...
           'gpuIDX.\n']);
    else
      wedgeType = 2;
      
      switch length(varargin) 
        case 7
          
          if (useV2)
            BH_templateSearch3d_2( varargin{2}, varargin{3},varargin{4}, ...
                              varargin{5}, varargin{6},wedgeType,varargin{7});
          else
            BH_templateSearch3d( varargin{2}, varargin{3},varargin{4}, ...
                              varargin{5}, varargin{6},wedgeType,varargin{7});
          end
 
        case 8
          if (useV2)    
       
            BH_templateSearch3d_2( varargin{2}, varargin{3},varargin{4}, ...
                                   varargin{5}, varargin{6},wedgeType, ...
                                   varargin{7},varargin{8});
          else
            BH_templateSearch3d( varargin{2}, varargin{3},varargin{4}, ...
                                 varargin{5}, varargin{6},wedgeType, ...
                                 varargin{7},varargin{8});
          end

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
  
  %%%%%%%% GLOBALS  %%%%%%%%%%%%%%%%%%%%

  % These variables are to maintain some flexibility for parameters that have
  % an ill-defined dependence on experimental factors. Preferably only until
  % they can be resolved.

  global bh_global_window_cutoff;
  try 
    bh_global_window_cutoff = pBH.('windowCutoff');
  catch
    bh_global_window_cutoff = -2;
  end
  
  % These are for making shape based masks. I think the problem is likely
  % dependent on the current resolution of the sub-tomogram, and that a
  % single set of values will not work for everything. 

  % Note that these must also be declared in the relevant functions
  
  %%%%%%% BH_mask3d.m %%%%%%%
  global bh_global_binary_mask_low_pass; 
  global bh_global_binary_mask_threshold;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%% BH_pcaPub.m %%%%%%%
  global bh_global_binary_pcaMask_threshold;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%% Anything that reads geometry. This way if size changes, its okay.
  %%%%%%% Needed only for tracking multiple copies of a single particle.
  %%%%%%% (Currently not used)
  global bh_global_nCol 
  bh_global_nCol = 26;
  %%%%%%%

  %%%%%% BH_mask3d - affects then the FSC calc
  global bh_global_vol_est_scaling;
  global bh_global_MTF;
  
  %%%%% 
  global bh_global_fast_scratch_disk;
  global bh_global_ram_disk;
  
  %%%%%%% BH_ctfCorrect_3d
  %%%%%%% Wiener filter and cut off past this point
  global bh_global_turn_on_phase_plate;

  try
    bh_global_turn_on_phase_plate = pBH.('phakePhasePlate');
  catch
    bh_global_turn_on_phase_plate = 0;
  end
  
  %%%%%%% BH_ctf_estimate, updateFFT
  %%%% Can't pad K3 images enough to avoid ghosting until mexInterp is
  %%%% ready
  global bh_global_do_2d_fourier_interp;
  try
    bh_global_do_2d_fourier_interp = pBH.('useFourierInterp');
  catch
    bh_global_do_2d_fourier_interp = 0;
  end
  
  global bh_global_save_tomoCPR_diagnostics;
  try
    bh_global_save_tomoCPR_diagnostics = pBH.('tomoCprDiagnostics');
  catch
    bh_global_save_tomoCPR_diagnostics = 0;
  end
  
  
%%%%%%%%%%%%%%

  %%%%% For profiling
  global bh_global_do_profile;
  try 
    bh_global_do_profile = pBH.('doProfile');
  catch
    bh_global_do_profile = false;
  end
    

  try
    bh_global_fast_scratch_disk  = pBH.('fastScratchDisk');
  catch
    bh_global_fast_scratch_disk='';
  end
  
  if ~isempty(bh_global_fast_scratch_disk)
    testFileName = sprintf('%s/thisEmCDiskCheck123456.txt',bh_global_fast_scratch_disk);
    [writeError] = system(sprintf('echo a > %s',testFileName));
    if (writeError)
      fprintf('\nRan into an error trying to write to the fastScatchDisk %s\n',bh_global_fast_scratch_disk);
    else
      system(sprintf('rm %s',testFileName));
      fprintf('Found and using your fastScratchDisk\n');
    end
  end
    
  try
    bh_global_ram_disk = pBH.('ramDisk');
  catch
    bh_global_ram_disk = '/dev/shm';
  end
 
  testFileName = sprintf('%s/thisEmCDiskCheck123456.txt',bh_global_ram_disk);
  [writeError] = system(sprintf('echo a > %s',testFileName));  
  if (writeError)
    fprintf('\nRan into an error trying to write to the fastScatchDisk %s\n',bh_global_ram_disk);
    bh_global_ram_disk = '';
  else
    fprintf('Found and using your ramDisk\n');
    system(sprintf('rm %s',testFileName));
  end

  
  try   
    bh_global_binary_mask_low_pass = pBH.('setMaskLowPass');
  catch
    % These seem to be okay for higher-resolution data (EMPIAR ribo sets)
    bh_global_binary_mask_low_pass = 14;
  end
  
  try
    bh_global_binary_mask_threshold = pBH.('setMaskThreshold');
  catch
    bh_global_binary_mask_threshold = 2.5;
  end
  
  try
    bh_global_binary_pcaMask_threshold = pBH.('setPcaMaskThreshold');
  catch
    bh_global_binary_pcaMask_threshold = 0.5;
  end

  global bh_global_kFactorScaling;
  try
    bh_global_kFactorScaling = pBH.('kFactorScaling');
  catch
    bh_global_kFactorScaling = 1.0;
  end
  
  global bh_global_tomoCPR_random_subset;
  try
    bh_global_tomoCPR_random_subset = pBH.('tomoCPR_randomSubset');
  catch
    bh_global_tomoCPR_random_subset = 500;
  end
  
  try
    bh_global_vol_est_scaling = pBH.('setParticleVolumeScaling');
  catch
    % The low pass version of the map used for the estimate overestimates
    % the molecular volume at the hydration radius of the underlying atoms.
    % This flag will override the value I've calculated which depends on
    % the masking resolution. TODO when the map resolution is lower than
    % the masking resolution, this will again underestimate the scaling,
    % artificialy *de*pressing the FSC. Left to zero this is calculated in
    % mask_3d
    bh_global_vol_est_scaling = 0.0;
  end
  
  try
    % 0 - off, 2 original (matches closely measured MTF), 1 stronger
    % Anthing else, float, iX = scalar, dX = cap val e.g. 
    % opiton 1 100.04 and 2 (default) is 25.06
    bh_global_MTF = pBH.('mtfVal');
  catch
    bh_global_MTF = 2;
  end

  global bh_global_print_shifts_in_particle_basis;
  try 
    bh_global_print_shifts_in_particle_basis = pBH.('printShiftsInParticleBasis');
  catch
    bh_global_print_shifts_in_particle_basis = true;
  end
  
  global bh_global_zero_lag_score;
  try 
    bh_global_zero_lag_score = pBH.('useZeroLagScore');
  catch
    bh_global_zero_lag_score = false;
  end
  
  global bh_global_ML_compressByFactor;
  global bh_global_ML_angleTolerance;
  try
    bh_global_ML_compressByFactor = pBH.('ML_compressByFactor')
  catch
    bh_global_ML_compressByFactor = 1.25
  end
  try
    bh_global_ML_angleTolerance = pBH.('ML_angleTolerance')
  catch
    bh_global_ML_angleTolerance = 5
  end
  
  fprintf('nExpGlobals %2.2f maskLP, %2.2f maskThr, %2.2f pcaMaskThr\n', ...
          bh_global_binary_mask_low_pass, ...
          bh_global_binary_mask_threshold, ...
          bh_global_binary_pcaMask_threshold);
catch
  error('error parsing parameter file %s\n', paramTest)
end



end

function [] = select_gpus(nGPUs_wanted,nGPUs_visible,cmdIN)

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


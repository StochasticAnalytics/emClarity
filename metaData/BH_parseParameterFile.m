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
last_parsed_parameter = 'none';
% Check that all paramters are name: value pairs
stringValues = {'subTomoMeta'; ...
  'Ali_mType';'Cls_mType';'Cls_mType';'Raw_mType';'Fsc_mType'; ...
  'Pca_distMeasure';'Kms_mType';'flgPrecision';'Tmp_xcfScale';...
  'fastScratchDisk';'Tmp_eraseMaskType';'startingDirection';'Peak_mType';'symmetry'};
for i = 1:size(p2,1)
  pNameVal = strsplit(p2{i,1},'=');
  if length(pNameVal) == 1
    fprintf("Last successfully parsed parameter: %s\n", string(last_parsed_parameter));
    error('Could not split Name=Value pair for\n\t %s',char(pNameVal))
  elseif length(pNameVal) > 2
    fprintf("Last successfully parsed parameter: %s\n", string(last_parsed_parameter));
    error('To many colons in\n\t %s',char(pNameVal))
  else   
    if any(strcmp(stringValues, pNameVal{1}))
      emc.(pNameVal{1}) = pNameVal{2};
    else
      emc.(pNameVal{1}) = EMC_str2double(pNameVal{2});
    end
    last_parsed_parameter = pNameVal{1};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Asserts on required parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(emc, 'nGPUs')
  EMC_assert_numeric(emc.nGPUs, 1, [1, 1000]);
else
  error('nGPUs is a required parameter');
end

if isfield(emc, 'nCpuCores')
  EMC_assert_numeric(emc.nCpuCores, 1, [1, 1000]);
else
  error('nCpuCores is a required parameter');
end

symmetry_has_been_checked = false;
if ~isfield(emc, 'symmetry')
  %TODO asserts on allowed values for symmetry paraemeter
  error('You must now specify a symmetry=X parameter, where symmetry E (C1,C2..CX,O,I)');
end
symmetry_has_been_checked = true;

if isfield(emc, 'PIXEL_SIZE')
  EMC_assert_numeric(emc.PIXEL_SIZE, 1, [0, 100e-10]);
  emc.pixel_size_si = emc.PIXEL_SIZE;
  emc.pixel_size_angstroms = emc.PIXEL_SIZE.*10^10;
else
  error('PIXEL_SIZE is a required parameter');
end

if isfield(emc, 'Cs')
  EMC_assert_numeric(emc.Cs, 1, [0, 10e-3]);
else
  error('Cs is a required parameter');
end

if isfield(emc, 'VOLTAGE')
  EMC_assert_numeric(emc.VOLTAGE, 1, [20e3, 1000e3]);
else
  error('VOLTAGE is a required parameter');
end

if isfield(emc, 'AMPCONT')
  EMC_assert_numeric(emc.AMPCONT, 1, [0.0, 1.0]);
  if emc.Cs == 0
    emc.Cs = 1e-10;
  end
  
else
  error('AMPCONT is a required parameter');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now check for optional parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Early development parameter, used to store more than one orientation during template matching
% and use for further refinement.
if isfield(emc, 'nPeaks')
    EMC_assert_numeric(emc.nPeaks, 1);
else
  emc.('nPeaks') = 1;
end

% Used when cutting out subtomos for further processing. Adds extra padding to anticipate shifts etc.
% This has not been well tested

% When used in average3d, this value is stored in the subTomoMeta.
if isfield(emc, 'CUTPADDING')
  EMC_assert_numeric(emc.CUTPADDING, 1);
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

if isfield(emc, 'flgCones')
  EMC_assert_boolean(emc.flgCones)
else
  emc.('Fsc_bfactor') = false;
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

if isfield(emc,'eucentric_fit')
  EMC_assert_boolean(emc.eucentric_fit)
else
  emc.eucentric_fit = false;
end

if isfield(emc,'eucentric_minTilt')
  EMC_assert_numeric(emc.eucentric_maxTilt)
else
  emc.eucentric_maxTilt = 50.0;
end


% TODO: these should maybe be two different orthogonal parameters
  % if > 1 keep this many subtomos
  % if < 1 keep this fraction
emc = EMC_assert_deprecated_substitution(emc, 0.0, 'ccc_cutoff', 'flgCCCcutoff');
EMC_assert_numeric(emc.ccc_cutoff,1)


% TOOD: DOC
emc = EMC_assert_deprecated_substitution(emc, false, 'projectVolumes', 'flgProjectVolumes');
EMC_assert_boolean(emc.projectVolumes);

% Whether the cycle is expected to be used for classification or alignment.
% Eventually, the distinction should not matter.
emc = EMC_assert_deprecated_substitution(emc, false, 'classification', 'flgClassify');
EMC_assert_boolean(emc.classification);


emc = EMC_assert_deprecated_substitution(emc, 0, 'multi_reference_alignment', 'flgMultiRefAlignment');
EMC_assert_numeric(emc.multi_reference_alignment, 1, [0, 2]);

% Zero padding of the volumes before alignment/other FFT ops
emc = EMC_assert_deprecated_substitution(emc, 1.5, 'scale_calc_size', 'scaleCalcSize');
EMC_assert_numeric(emc.scale_calc_size, 1, [1.0, 2.0]);

emc = EMC_assert_deprecated_substitution(emc, false, 'limit_to_one_core', 'flgLimitToOneProcess');
EMC_assert_boolean(emc.limit_to_one_core);

if (emc.limit_to_one_core)
  emc.nCpuCores = 1;
end

if isfield(emc, 'force_no_symmetry')
  EMC_assert_boolean(emc.force_no_symmetry)
  if (~symmetry_has_been_checked)
    error('force_no_symmetry must be after symmetry check');
  end
  % Warning must be after symmetry check
  if (emc.force_no_symmetry)
    emc.symmetry='C1';
  end
else
  emc.force_no_symmetry = true;
end

if isfield(emc, 'Pca_constrain_symmetry')
  EMC_assert_boolean(emc.Pca_constrain_symmetry)
else
  emc.Pca_constrain_symmetry = false;
end

emc = EMC_assert_deprecated_substitution(emc, false, 'fsc_with_chimera', 'fscWithChimera');
EMC_assert_boolean(emc.fsc_with_chimera);

emc = EMC_assert_deprecated_substitution(emc, 0.1, 'minimum_particle_for_fsc_weighting', 'minimumparticleVolume');
EMC_assert_numeric(emc.minimum_particle_for_fsc_weighting, 1, [0.01, 1.0]);

emc = EMC_assert_deprecated_substitution(emc, 1.0, 'fsc_shape_mask', 'flgFscShapeMask');
EMC_assert_numeric(emc.fsc_shape_mask, 1, [0.0, 2.0]);

if isfield(emc, 'shape_mask_lowpass')
  EMC_assert_numeric(emc.shape_mask_lowpass, 1, [10, 100]);
else
  emc.shape_mask_lowpass = 14;
end

if isfield(emc, 'shape_mask_threshold')
  EMC_assert_numeric(emc.shape_mask_threshold, 1, [0.1, 10.0]);
else
  emc.shape_mask_threshold = 2.4;
end

if isfield(emc, 'shape_mask_test')
  EMC_assert_boolean(emc.shape_mask_test);
else
  emc.shape_mask_test = false;
end

emc = EMC_assert_deprecated_substitution(emc, 22.0, 'pca_scale_spaces', 'pcaScaleSpace');
EMC_assert_numeric(emc.pca_scale_spaces);
emc.('n_scale_spaces') = numel(emc.pca_scale_spaces);

if isfield(emc, 'Pca_maxEigs')
  EMC_assert_numeric(emc.Pca_maxEigs, 1, [1, 1000]);
else
  emc.Pca_maxEigs = 36;
end

if isfield(emc, 'Pca_randSubset')
  EMC_assert_numeric(emc.Pca_randSubset, 1);
else
  emc.Pca_randSubset = 0;
end

if ~isfield(emc, 'Pca_clusters');
  error('Pca_clusters is a required parameter');
end

% Allowed values are validated inside BH_clusterPub.m
if ~isfield(emc, 'Pca_distMeasure');
  emc.distance_metric = 'sqeuclidean';
end

if isfield(emc, 'Pca_nReplicates');
  EMC_assert_numeric(emc.Pca_nReplicates, 1, [100, 1000]);
else
  emc.n_replicates = 256;
end

if isfield(emc, 'Pca_refineKmeans')
  EMC_assert_boolean(emc.Pca_refineKmeans);
else
  emc.Pca_refineKmeans = false;
end

if isfield(emc, 'Pca_flattenEigs')
  EMC_assert_boolean(emc.Pca_flattenEigs);
else
  emc.Pca_flattenEigs = true;
end


if isfield(emc, 'Pca_som_coverSteps')
  EMC_assert_numeric(emc.Pca_som_coverSteps, 1, [1, 1000]);
else
  emc.Pca_som_coverSteps = 100;
end

if isfield(emc, 'Pca_som_initNeighbor')
  EMC_assert_numeric(emc.Pca_som_initNeighbor, 1, [1, 32]);
else
  emc.Pca_som_initNeighbor = 3;
end

if ~isfield(emc, 'Pca_som_topologyFcn')
  % TODO: assert on allowed values
  emc.Pca_som_topologyFcn = 'hextop';
end

if isfield(emc, 'spike_prior')
  EMC_assert_boolean(emc.spike_prior);
else
  emc.spike_prior = false;
end


emc = EMC_assert_deprecated_substitution(emc, false, 'update_class_by_ccc', 'updateClassByBestReferenceScore');
EMC_assert_boolean(emc.update_class_by_ccc);
if (~emc.multi_reference_alignment)
  % update by ccc only makes sense for multi reference alignment
  emc.update_class_by_ccc = false;
end

emc = EMC_assert_deprecated_substitution(emc, true, 'move_reference_by_com', 'flgCenterRefCOM');
EMC_assert_boolean(emc.move_reference_by_com);

if isfield(emc, 'use_new_grid_search')
  EMC_assert_boolean(emc.use_new_grid_search);
else
  emc.use_new_grid_search = true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%% tomoCPR params, mostly experimental

emc = EMC_assert_deprecated_substitution(emc, false, 'save_mapback_classes', 'flgColorMap');
EMC_assert_boolean(emc.save_mapback_classes);


% These seemed to be necessary at some point to translate between IMOD and emClarity 
% coordinate systems, but the should probably be looked at again. TODO:
if isfield(emc, 'flgPreShift')
  EMC_assert_numeric(emc.flgPreShift, 3);
else
  emc.flgPreShift = [-0.5,-0.5,0.5];
end


% These seemed to be necessary at some point to translate between IMOD and emClarity 
% coordinate systems, but the should probably be looked at again. TODO:
if isfield(emc, 'flgPostShift')
  EMC_assert_numeric(emc.flgPostShift, 2);
else
  emc.flgPostShift = [-0.5,-0.5];
end

if isfield(emc, 'prjVectorShift')
  EMC_assert_numeric(emc.prjVectorShift, 3);
else
  emc.prjVectorShift = [0.5,0.5,1.0]';
end

if isfield(emc,'pixelShift')
  EMC_assert_numeric(emc.pixelShift, 1);
else
  emc.pixelShift = -1;
end

if isfield(emc, 'pixelMultiplier')
  EMC_assert_numeric(emc.pixelMultiplier, 1);
else
  emc.pixelMultiplier = 1;
end

if isfield(emc, 'tomoCPR_random_subset')
  EMC_assert_numeric(emc.tomoCPR_random_subset, 1);
else
  emc.tomoCPR_random_subset = -1;
end

% I think this has been removed
if isfield(emc, 'probabilityPeakiness')
  EMC_assert_numeric(emc.probabilityPeakiness, 1);
else
  emc.probabilityPeakiness = 0;
end

if isfield(emc, 'whitenProjections')
  EMC_assert_numeric(emc.whitenProjections, 1);
else
  emc.whitenProjections = 0;
end

if isfield(emc, 'rot_option_global')
  EMC_assert_numeric(emc.rot_option_global, 1);
else
  emc.rot_option_global = 1;
end


if isfield(emc, 'rot_option_local')
  EMC_assert_numeric(emc.rot_option_local, 1);
else
  emc.rot_option_local = 1;
end

if isfield(emc, 'rot_default_grouping_global')
  EMC_assert_numeric(emc.rot_default_grouping_global, 1);
else
  emc.rot_default_grouping_global = 3;
end

if isfield(emc, 'rot_default_grouping_local')
  EMC_assert_numeric(emc.rot_default_grouping_local, 1);
else
  emc.rot_default_grouping_local = 3;
end

if isfield(emc, 'mag_option_global')
  EMC_assert_numeric(emc.mag_option_global, 1);
else
  emc.mag_option_global = 1;
end

if isfield(emc, 'mag_option_local')
  EMC_assert_numeric(emc.mag_option_local, 1);
else
  emc.mag_option_local = 1;
end

if isfield(emc, 'mag_default_grouping_global')
  EMC_assert_numeric(emc.mag_default_grouping_global, 1);
else
  emc.mag_default_grouping_global = 5;
end

if isfield(emc, 'mag_default_grouping_local')
  EMC_assert_numeric(emc.mag_default_grouping_local, 1);
else
  emc.mag_default_grouping_local = 5;
end

if isfield(emc, 'tilt_option_global')
  EMC_assert_numeric(emc.tilt_option_global, 1);
else
  emc.tilt_option_global = 5;
end

if isfield(emc, 'tilt_option_local')
  EMC_assert_numeric(emc.tilt_option_local, 1);
else
  emc.tilt_option_local = 5;
end

if isfield(emc, 'tilt_default_grouping_global')
  EMC_assert_numeric(emc.tilt_default_grouping_global, 1);
else
  emc.tilt_default_grouping_global = 5;
end

if isfield(emc, 'tilt_default_grouping_local')
  EMC_assert_numeric(emc.tilt_default_grouping_local, 1);
else
  emc.tilt_default_grouping_local = 5;
end


if isfield(emc, 'peak_mask_fraction')
  EMC_assert_numeric(emc.flgPeakMask,1);
else
  emc.peak_mask_fraction = 0.4;
end

if isfield(emc, 'min_overlap')
  EMC_assert_numeric(emc.min_overlap,1);
else
  emc.min_overlap = 0.5;
end


if isfield(emc, 'k_factor_scaling')
  EMC_assert_numeric(emc.k_factor_scaling,1);
else
  emc.k_factor_scaling = nan;
end

if isfield(emc, 'shift_z_to_to_centroid')
  EMC_assert_boolean(emc.shift_z_to_to_centroid);
else
  emc.shift_z_to_to_centroid = true;
end

emc = EMC_assert_deprecated_substitution(emc, 500e-9, 'tomo_cpr_defocus_range', 'tomoCprDefocusRange');
EMC_assert_numeric(emc.tomo_cpr_defocus_range, 1, [0.0, 10000e-9]);

emc = EMC_assert_deprecated_substitution(emc, 100e-9, 'tomo_cpr_defocus_step', 'tomoCprDefocusStep');
EMC_assert_numeric(emc.tomo_cpr_defocus_step, 1, [1.0e-9, 10000e-9]);

emc = EMC_assert_deprecated_substitution(emc, false, 'tomo_cpr_defocus_refine', 'calcCTF');
EMC_assert_boolean(emc.tomo_cpr_defocus_refine);

end


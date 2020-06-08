% This is a comment
% Inline comments will break the parser.


% String to name the structure that contains all of the metadata, projectName
subTomoMeta=emClarity_tutorial

fastScratchDisk=/scratch/himesb

% Number of GPUS
nGPUs=2
nCpuCores=12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    Mask parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The particle radius in x,y,z Angstrom, smallest value to contain particle. 
% For particles in a lattice, neighboring particles can be used in alignment
% by specifying a larger mask size, but this paramter must correspond to your
% target, a cetral hexamer of capsid proteins for example.
particleRadius=[160,160,160]
% Estimated particle mass in megaDa
particleMass=4

Ali_mType=sphere
Cls_mType=sphere

% For special cases where repeated motifs are present which might cause one
% subtomo to drift to a neighbor. This allows a larger alignment mask to be used
% for the rotational search (Ali_m...) but limits the translational peak search.
%Peak_mType=sphere
%Peak_mRadius=[100,100,100]

% mask radius and center - and center in Angstrom. Mask size is determined 
% large enough to contain delocalized signal, proper apodization, and to 
% avoid wraparound error in cross-correlation.
Ali_mRadius=[180,180,180]
Ali_mCenter=[0,0,0]
Cls_mRadius=[180,180,180]
Cls_mCenter=[0,0,0 ]


% Sampling rate
Ali_samplingRate=4
Cls_samplingRate=4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Tomo-constrained projection refinement parameters    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% I advise to avoid using this experimental feature for now.
tomoCprDefocusRefine=0
tomoCprDefocusRange=500e-9; 
tomoCprDefocusStep=20e-9;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    subTomogram alignment           %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Raw_className=0
% Second row specifies C1 symmetry
Raw_classes_odd=[0;1.*ones(2,1)]
Raw_classes_eve=[0;1.*ones(2,1)]

Raw_angleSearch=[15,5,0,0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Template matching parameters    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Tmp_samplingRate=6
Tmp_threshold=400
Tmp_angleSearch=[180,15,180,15,0]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Class reference   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




Cls_className=9
Cls_classes_odd=[1:9;1.*ones(1,9)]
Cls_classes_eve=[1:9;1.*ones(1,9)]

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    FSC Paramters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On/Off anisotropic SSNR calc
flgCones=1

% B-factor applied to weighted averages and refs. Should be < 20. Can be a vector
% where the 2:end positions generate independent maps at that sharpening 
% when "avg paramN.m N FinalAlignment is run.
Fsc_bfactor=10

% For very tightly packed subTomos set to 1 to avoid mixing halfsets
% form overlaping peripheral density.
fscGoldSplitOnTomos=0

% A very soft mask based on the particle envelope. Highly recommended to not 
% turn this off. If you suspect some density is being ommited, you could try 
% lowering the threshold. If you have extra "dust" outside your particle, it
% can be helpful to decrease the resolution used in the inital thresholding
% (shape_mask_lowpass) or to increase the threshold (or both).
% shape_mask_test will create the mask and save it to your FSC directory, which 
% is useful for testing these parameters. >$ emClarity fsc paramN.m N RawAlignment
flgFscShapeMask=1
shape_mask_lowpass=18
shape_mask_threshold=2.4
shape_mask_test=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Classification Paramters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On/Off classification. This must be on when "avg paramN.m N RawAlignment"
% is run at the begining of a cycle where classification is to be run.
flgClassify=0

% List of different cluster sizes to try, eg [3;4]
Pca_clusters=[2,3]

% Maximum number of eigenvalues/vectors to save
Pca_maxEigs=36

% Different resolution bands to run PCA on. Not all need to be used for subsequent
% clustering. (Angstrom)
pcaScaleSpace=[8,14,21];

% Random subset of particles used to reduce the burden of PCA
% This is ignored if flgPcaFull is true in the call to "pca"
Pca_randSubset=0

% Different ranges of coefficients to use in the clustering. At times, the 
% missing wedge can be a strong feature, such that ignoring the first few 
% eigen values can be usefule. [2:40 ; 6;40 ; 10:40]
% Each row must have the same number of entries, and there must be a row 
% for each scale space, even if it is all zeros.
Pca_coeffs=[7:11;7:11;7:11]



% The number of subtomos to process at once before pulling tempDataMatrix off 
% the gpu and into main memory.
PcaGpuPull=1200





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for CTF all si (meters, volts)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%   Microscope parameters     %%%%%%%%%%

% Of the data saved in fixed stacks - MUST match header
PIXEL_SIZE=2.17e-10 
% Currently any super-resolution data is cropped in Fourier Space after alignment
% allowing for finer sampling when interpolating the stacks, while then 
% filtering out noise due to aliasing.
SuperResolution=0
% Spherical abberation
Cs=2.7e-3 
% Accelerating voltage
VOLTAGE=300e3
% Percent amplitude contrast
AMPCONT=0.10

% search range - generally safe to test a wide range
defEstimate=3.5e-6
defWindow=2e-6
% The PS is considered from the lower resolution inflection point
% past the first zero to this cutoff resolution
defCutOff=7e-10

% Total dose in electron/A^2, assumed constant rate
CUM_e_DOSE=60
% Gold fiducial diameter
beadDiameter=10e-9



% Boolean: If a graduated dose scheme is used, where applied dose is ~ Thickness (1/cos(tiltAngle))
oneOverCosineDose=0
% Degrees: first tilt collected
startingAngle=30
% String: Direction after first tilt (pos or neg)
startingDirection=neg
% Boolean: Group size on either side of zero, e.g. 2 (with pos) would give, 0, 3, 6, -3, -6, 9, ...
doseSymmetricIncrement=0
% Float: Dose applied per tilt (at zero degrees)
doseAtMinTilt=1.46

% This is a comment
% Inline comments will break the parser.


% String to name the structure that contains all of the metadata, projectName
subTomoMeta=emClarity_tutorial

% Number of GPUS, either a list [1,3] or one int indicating a range 1:int
nGPUs=4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    Mask parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The particle radius in x,y,z Angstrom, smallest value to contain particle. 
% For particles in a lattice, neighboring particles can be used in alignment
% by specifying a larger mask size, but this paramter must correspond to your
% target, a cetral hexamer of capsid proteins for example.
particleRadius=[160,160,160]

Ali_mType=sphere
Cls_mType=sphere

% large enough to contain delocalized signal, proper apodization, and to 
% avoid wraparound error in cross-correlation.
Ali_mRadius=[180,180,180]
Ali_mCenter=[ 0,0,0 ]
Cls_mRadius=[180,180,180]
Cls_mCenter=[ 0,0,0 ]


% Sampling rate
Ali_samplingRate=3
Cls_samplingRate=3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Tomo-constrained projection refinement parameters    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tomoCprLowPass=18


tomoCprDefocusRefine=1
tomoCprDefocusRange=300e-9; 
tomoCprDefocusStep=15e-9;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    subTomogram alignment           %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Raw_className=0
Raw_classes_odd=[0;1.*ones(2,1)]
Raw_classes_eve=[0;1.*ones(2,1)]

Raw_angleSearch=[0,0,12,1.5]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Template matching parameters    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Tmp_samplingRate=5
Tmp_threshold=400
Tmp_angleSearch=[180,15,180,15,0]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Class reference and alignment    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NoA_className=0 
NoA_classes_odd=[0;1.*ones(1,1)]
NoA_classes_eve=[0;1.*ones(1,1)]

ref_AngleShift_odd=0
ref_TransShift_odd=[0,0,0]
ref_Ref_odd=[1;1]
ref_Sharpen=0

ref_AngleShift_eve=0
ref_TransShift_eve=[0,0,0]
ref_Ref_eve=[1;1]

Ref_angleSearch=[3,1.5,6,1;3,1.5,6,1]

Ref_className=4
Ref_classes_odd=[1:4;1.*ones(1,4)]
Ref_classes_eve=[1:4;1.*ones(1,4)]

Ref_references_odd=[1:6;1.*ones(1,6);ones(1,6)]
Ref_references_eve=[1:6;1.*ones(1,6);ones(1,6)]

Cls_className=9
Cls_classes_odd=[1:9;1.*ones(1,9)]
Cls_classes_eve=[1:9;1.*ones(1,9)]



Cls_angleSearch=[3,1,4,0.5]

% The out of plane inc must be a factor of the range, however, the in plane does 
% not need to be, and should be chosen to avoid symmetry related peaks.

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    FSC Paramters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% References to compare, sorted by order
Fsc_refVector_odd=[1]
Fsc_refVector_eve=[1]


Fsc_bfactor=40

% For very tightly packed subTomos set to 1 to avoid mixing halfsets
% form overlaping peripheral density.
fscGoldSplitOnTomos=0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Classification Paramters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of different cluster sizes to try, eg [3;4]
Pca_clusters=[2,3]

% Maximum number of eigenvalues/vectors to save
Pca_maxEigs=36
pcaScaleSpace=[8,14,21];
% Random subset of particles used to reduce the burden of SVD
% This is ignored if flgPcaFull is true in the call to BH_pca
Pca_randSubset=0
% Clustering parameters
Pca_flattenEigs=0;
Pca_relativeScale=[1,1,1]
% Different ranges of coefficients to use in the clustering. At times, the 
% missing wedge can be a strong feature, such that ignoring the first few 
% eigen values can be usefule. [2:40 ; 6;40 ; 10:40]
Pca_coeffs_odd=[7,0,0;14,15,0;1,19,20]
Pca_coeffs_eve=[7,0,0;14,15,0;1,19,20]


% distance measure for k-means. I've decided cosine is marginally better 
% for my data, feel free to test. sqeuclidean, cityblock, cosine, correlation.
Pca_distMeasure=sqeuclidean
Pca_nReplicates=128
Pca_nCores=24
% The number of subtomos to process at once before pulling tempDataMatrix off 
% the gpu and into main memory.
PcaGpuPull=1200
flgWedgeComparison=0
flgWMDs=3





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for CTF all si (meters, volts)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%   3d CTF paramters    %%%%%%%%%%

% a fraction of the average X,Y dimension to minimally construct to allow for 
% accurate r-weighting. 
flgSlabPad=0.25
% dimension perpendicular to the tilt-axis to pad to during ctf multiplication.
flgCtfStripWidth=1536
% Used to deterine slab height and strip width during ctf correction
defocusErrorEst=50e-9


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

% test astigmatism vals (run or not, meters, degrees)
flgAstigmatism=1
maxAstig=2000e-10
astigStep=150e-10
coarseAngStep=10
fineAngStep=0.5

% numWorkers only affects the interpolation of the stack (spline on cpu)

% The following can be adjusted but be careful

% Change in z height allowed for tile inclusion in global estimate
deltaZTolerance=50e-9
% Offset, use to determine tilt gradient.
zShift=0
% Tile size, smaller = smoother, but longer runs. Larger values needed to ensure
% high frequency information is included.
tileSize=320
% overlap as a fraction of tilesize (tileSize/tileOverlap)
tileOverlap=4
% Padded size to reduce aliasing effects, more important as defocus becomes higher.
% for defocus > 6um use 1536, 1024 prob okay for lower. On any decent card, going smaller
% doesn't really speed things up so don't bother
paddedSize=1536



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Advanced options    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flgPrecision=single
flgCones=2
flgClassify=0
flgSymmetrizeSubTomos=0
duplicateRadius=8
duplicateSampling=4
flgCCCcutoff=0.00
removeBottomPercent=0.0
flgMultiRefAlignment=0
flgRotAvgRef=0

experimentalOpts=[2,1,1,1.5,4,1,1,3,1,0,1]
fastScratchDisk=

rmsScale=1.0

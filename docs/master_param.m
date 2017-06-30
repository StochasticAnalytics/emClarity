% This is a comment
% Inline comments will break the parser.

%%%% TODO
% pixel size is specified in three places --> one
% Simplify mVals/SamplingRate


% String to name the structure that contains all of the metadata, projectName
subTomoMeta=rln_tutorial_1


% Number of GPUS, either a list [1,3] or one int indicating a range 1:int
nGPUs=1
% 1 use as much memory as possible | 2 conserve memory. Always try to run in 1
% but often at full sampling 2 is necessary.
avgMemory=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    Mask parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoA_mType=sphere
Pca_mType=sphere
Cls_mType=sphere
Raw_mType=sphere
Fsc_mType=sphere
Kms_mType=sphere

Fsc_shapeMask=0
Pca_shapeMask=0
Raw_shapeMask=0

% mask size, radius, center - note here mask size is ignored, and the window size
% is used insead. Note that each is a 3 x 3 matrix
NoA_mVals=[156,156,156; 84,84,84; 0,0,0]
Pca_mVals=[156,156,156; 84,84,84; 0,0,0]
Cls_mVals=[156,156,156; 84,84,84; 0,0,0]
Raw_mVals=[156,156,156; 84,84,84; 0,0,0]
Fsc_mVals=[156,156,156; 84,84,84; 0,0,0]
Kms_mVals=[156,156,156; 84,84,84; 0,0,0]


% Sampling rate
% Note that all pixel values should be entered at full binning 
NoA_samplingRate=4
Pca_samplingRate=4
Cls_samplingRate=4
Raw_samplingRate=4
Kms_samplingRate=4
Fsc_samplingRate=4
Ref_samplingRate=4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Advanced options    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


scaleCalcSize=2
flgReWeight=1
flgQualityWeight=1
flgPrecision=single
interpolationOrder=1
flgCones=0
flgClassify=0
flgMultiRefAlignment=0
flgRotAvgRef=0
flgCenterRefCOM=1
duplicateRadius=8
duplicateSampling=3
flgCCCcutoff=0.00
removeBottomPercent=0.0


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

Cls_peakSearch=[30,30,24]
Cls_centerOfMass=[2,2,2]
Cls_angleSearch=[3,1,4,0.5]

% The out of plane inc must be a factor of the range, however, the in plane does 
% not need to be, and should be chosen to avoid symmetry related peaks.

 
Raw_className=0
Raw_classes_odd=[0;1.*ones(2,1)]
Raw_classes_eve=[0;1.*ones(2,1)]

Raw_peakSearch=[64,64,64]
Raw_centerOfMass=[3,3,3]
Raw_angleSearch=[0,0,12,3,0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    FSC Paramters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fsc_withChimera=1
% References to compare, sorted by order
Fsc_refVector_odd=[1]
Fsc_refVector_eve=[1]
% Alignment parameters to bring the two half sets into register
Fsc_peakSearch=[84,84,84]
Fsc_centerOfMass=[3,3,3]
Fsc_angleSearch=[3,0.5,3,0.3]


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
%%%%%%%%%%%%%%%%%%%    Template matching parameters    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Tmp_samplingRate=6
Tmp_threshold=400
Tmp_targetSize=[512,512,512]

% lattice param is ~ 120A for trimer, 69A for unitcell
Tmp_latticeRadius=[160,160,160]
Tmp_angleSearch=[180,12,180,12,0]
Tmp_xcfScale=xcf

xfcSize=[256,256,256]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Tomo-particle polishing parameters    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mapBackIter=0
mapBackLowPass=12
mapBackRePrjSize=128



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for CTF all si (meters, volts)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%   3d CTF paramters    %%%%%%%%%%

% a fraction of the average X,Y dimension to minimally construct to allow for 
% accurate r-weighting. 
flgSlabPad=0.25
% dimension perpendicular to the tilt-axis to pad to during ctf multiplication.
flgCtfStripWidth=2048
% Used to deterine slab height and strip width during ctf correction
defocusErrorEst=30e-9


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
AMPCONT=0.085

% search range - generally safe to test a wide range
defEstimate=4e-6
defWindow=2.5e-6
% The PS is considered from the lower resolution inflection point
% past the first zero to this cutoff resolution
defCutOff=8e-10

% Total dose in electron/A^2, assumed constant rate
CUM_e_DOSE=100
% Gold fiducial diameter
beadDiameter=10e-9

% test astigmatism vals (run or not, meters, degrees)
flgAstigmatism=1
maxAstig=2000e-10
astigStep=150e-10
coarseAngStep=10
fineAngStep=0.5

% numWorkers only affects the interpolation of the stack (spline on cpu)
numWorkers=7

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

%%%%%%%%%%%%%%%%%%    Developmental Parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rln_mVals=[160,160,120; 48,48,40; 0,0,0]
Rln_samplingRate=1
Rln_tomogramCp=0
% 0 - nothing,1 - mv, 2 - cp , 3 - ln symbolic


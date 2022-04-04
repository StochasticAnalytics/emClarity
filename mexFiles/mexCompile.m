function [] = mexCompile(varargin)

fprintf("\n\nCompile here\n\n");
mexPATH = '/groups/himesb/git/emClarity/mexFiles/';
CUDA_LIB = '-L/groups/cryoadmin/software/CUDA-TOOLKIT/cuda_11.6.0/lib64';   ... % NOTE if you leave a space at the end of this string, MATLAB does not parse the option correctly (which wouldn't matter in a normal compile line!)
setenv('MW_NVCC_PATH','/groups/cryoadmin/software/CUDA-TOOLKIT/cuda_11.6.0/bin');
setenv('MW_ALLOW_ANY_CUDA','yes');
% getenv('CUDA_HOME')

system(sprintf('mkdir -p %s', mexPATH));
% For now just included everything in total.
inc = {'rotation_matrix.cpp','ctf.cu'};
for i = 1:length(inc)
  inc{i} = sprintf('%sutils/%s',mexPATH,inc{i});
end
% mexFILE = 'mexFFT' 'mexXform2d' 'mexCTF',;
if nargin > 0
  mexFILE = varargin;
else
 mexFILE = {'mexCTF','mexFFT','mexXform3d','mexSF3D'};
%  mexFILE = {'mexXform3d'};
end
% --extra-device-vectorization 
% --restrict
% --warn-on-double-precision-use
% --warn-on-spills
% -Wno-deprecated-gpu-targets
% add '-v' to troubleshood mexcuda compilatoin
mexcuda_opts = { ...
CUDA_LIB ...
'-lcublas'          ...            % Link to cuBLAS
'-lmwlapack'        ...            % Link to LAPACK
'-lcufft'           ...            % Link to cuFFT
['NVCCFLAGS=  --use_fast_math --default-stream per-thread -m64 --extra-device-vectorization  '...
 '--gpu-architecture=compute_80 ' ...
 '--restrict -Xptxas --warn-on-spills ' ...
 '-gencode=arch=compute_70,code=sm_70 ' ...
 '-gencode=arch=compute_75,code=sm_75 ' ...
 '-gencode=arch=compute_80,code=sm_80 ' ...
 '-gencode=arch=compute_86,code=sm_86 ' ...
 '-gencode=arch=compute_87,code=sm_87 '] ...% the optimizations are default anyway when I checked 

};

% '-L/usr/local/cuda-9.1/lib64'   ...    % Location of CUDA libraries
% '-L/usr/local/cuda-9.1/nvvm/lib64'};
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/lib64'   ...    % Location of CUDA libraries
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/nvvm/lib64'};


system('nvcc --version');
for i =1: length(mexFILE)
  mexcuda( mexcuda_opts{:}, sprintf('%s%s.cu',mexPATH,mexFILE{i}), inc{1}, inc{2});
  system('pwd');
  system(sprintf('mv %s.mexa64 %s/compiled',mexFILE{i}, mexPATH));
end
end

%nvcc -ptx --library-path /groups/grigorieff/home/himesb/thirdParty/cuda-9.2/nvvm/lib64




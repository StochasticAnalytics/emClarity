function [] = mexCompile(varargin)

fprintf("\n\nCompile here\n\n");
mexPATH = '/scratch/salina/git/emClarity/mexFiles/';
CUDA_LIB = '-L/usr/local/cuda/lib64';   ... % NOTE if you leave a space at the end of this string, MATLAB does not parse the option correctly (which wouldn't matter in a normal compile line!)
getenv('MW_NVCC_PATH')
getenv('CUDA_HOME')

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
mexcuda_opts = { ...
CUDA_LIB ...
'-lcublas'          ...            % Link to cuBLAS
'-lmwlapack'        ...            % Link to LAPACK
'-lcufft'           ...            % Link to cuFFT
['NVCCFLAGS=  --use_fast_math --default-stream per-thread -m64 '...
 '--gpu-architecture=compute_86 ' ...
 '--restrict -Xptxas --warn-on-spills ' ...
 '-gencode=arch=compute_60,code=sm_70 ' ...
 '-gencode=arch=compute_61,code=sm_75 ' ...
 '-gencode=arch=compute_70,code=sm_80 ' ...
 '-gencode=arch=compute_75,code=sm_86 '] ...% the optimizations are default anyway when I checked 

};




for i =1: length(mexFILE)
  mexcuda( mexcuda_opts{:}, sprintf('%s%s.cu',mexPATH,mexFILE{i}), inc{1}, inc{2});

  system(sprintf('mv %s.mexa64 %s/compiled',mexFILE{i}, mexPATH));
end
end





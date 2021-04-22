function [] = mexCompile(varargin)

fprintf("\n\nCompile here\n\n");
mexPATH = '~/emClarity/mexFiles/';
%setenv('nvcc','/admin/software/cuda-10.2/bin/nvcc');

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
'-L/misc/local/cuda-10.0/lib64 '   ... 
'-L/misc/local/cuda-10.0/nvvm/lib64 ' ...% Location of CUDA libraries
'-lcublas'          ...            % Link to cuBLAS
'-lmwlapack'        ...            % Link to LAPACK
'-lcufft'           ...            % Link to cuFFT
['NVCCFLAGS=  --use_fast_math --default-stream per-thread -m64 '...
 '--gpu-architecture=compute_75 ' ...
 '--restrict -Xptxas --warn-on-spills ' ...
 '-gencode=arch=compute_60,code=sm_60 ' ...
 '-gencode=arch=compute_61,code=sm_61 ' ...
 '-gencode=arch=compute_70,code=sm_70 ' ...
 '-gencode=arch=compute_75,code=sm_75 '] ...% the optimizations are default anyway when I checked 
};

% '-L/usr/local/cuda-9.1/lib64'   ...    % Location of CUDA libraries
% '-L/usr/local/cuda-9.1/nvvm/lib64'};
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/lib64'   ...    % Location of CUDA libraries
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/nvvm/lib64'};



for i =1: length(mexFILE)
  mexcuda( mexcuda_opts{:}, sprintf('%s%s.cu',mexPATH,mexFILE{i}), inc{1}, inc{2});

  system(sprintf('mv %s.mexa64 %s/compiled',mexFILE{i}, mexPATH));
end
end

%nvcc -ptx --library-path /groups/grigorieff/home/himesb/thirdParty/cuda-9.2/nvvm/lib64




function [] = mexCompile(varargin)

fprintf("\n\nCompile here\n\n");
mexPATH = '~/work/emClarity/mexFiles/';

system(sprintf('mkdir -p %s', mexPATH));
% For now just included everything in total.
inc = {'rotation_matrix.cpp','ctf.cu'};
for i = 1:length(inc)
  inc{i} = sprintf('%sutils/%s',mexPATH,inc{i});
end
% mexFILE = 'mexFFT';
if nargin > 0
  mexFILE = varargin;
else
 mexFILE = {'mexCTF','mexFFT','mexXform3d','mexXform2d'};
%  mexFILE = {'mexXform3d'};
end

mexcuda_opts = { ...
'-lcublas'          ...            % Link to cuBLAS
'-lmwlapack'        ...            % Link to LAPACK
'-lcufft'           ...            % Link to cuFFT
['NVCCFLAGS= -Wno-deprecated-gpu-targets --use_fast_math --default-stream per-thread -m64 --gpu-architecture=compute_70 --gpu-code=sm_70,sm_75'] ...% the optimizations are default anyway when I checked 
'-L/groups/grigorieff/home/himesb/thirdParty/cuda-10.1.2/cuda-toolkit/lib64'   ...    % Location of CUDA libraries
'-L/groups/grigorieff/home/himesb/thirdParty/cuda-10.1.2/cuda-toolkit/lib64'};

% '-L/usr/local/cuda-9.1/lib64'   ...    % Location of CUDA libraries
% '-L/usr/local/cuda-9.1/nvvm/lib64'};
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/lib64'   ...    % Location of CUDA libraries
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/nvvm/lib64'};



for i =1: length(mexFILE)
  
  mexcuda(mexcuda_opts{:}, sprintf('%s%s.cu',mexPATH,mexFILE{i}), inc{1}, inc{2});

  system(sprintf('mv %s.mexa64 %s/compiled',mexFILE{i}, mexPATH));
end
end

%nvcc -ptx --library-path /groups/grigorieff/home/himesb/thirdParty/cuda-9.2/nvvm/lib64




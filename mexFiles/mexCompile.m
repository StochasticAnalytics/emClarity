function [] = mexCompile(varargin)

fprintf("\n\nCompile here\n\n");
mexPATH = '/sa_shared/git/emClarity/mexFiles';
CUDA_LIB = '-L/usr/local/cuda/lib64';   ... % NOTE if you leave a space at the end of this string, MATLAB does not parse the option correctly (which wouldn't matter in a normal compile line!)
getenv('MW_NVCC_PATH')
getenv('CUDA_HOME')

system(sprintf('mkdir -p %s', mexPATH));
% For now just included everything in total.
inc = {'rotation_matrix.cpp','ctf.cu'};
for i = 1:length(inc)
  inc{i} = sprintf('%s/utils/%s',mexPATH,inc{i});
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
% '-lcublas_static'          ...            % Link to cuBLAS

mexcuda_opts = { ...
CUDA_LIB ...
'-lcuda' ...
'-lmwlapack'        ...            % Link to LAPACK
'-lcufft' ...
'-lculibos' ...
'-lcudart' ...
'-ldl' ...
'-lrt'           ...            % Link to cuFFT
['NVCCFLAGS=  --use_fast_math --default-stream per-thread -m64  --extra-device-vectorization --expt-relaxed-constexpr -t8 '...
 '--gpu-architecture=compute_86 ' ...
 '--restrict -Xptxas --warn-on-spills ' ...
 '-gencode=arch=compute_70,code=sm_70 ' ...
 '-gencode=arch=compute_80,code=sm_80 ' ...
 '-gencode=arch=compute_75,code=sm_75 ' ...
 '-gencode=arch=compute_86,code=sm_86 ' ...
 '-gencode=arch=compute_89,code=sm_89 '] ...% the optimizations are default anyway when I checked 

};

disp(mexcuda_opts);
if isfile(sprintf('%s/compiled',mexPATH))
  system(sprintf('rm -rf %s/compiled',mexPATH));
  system(sprintf('mkdir -p %s/compiled',mexPATH));
end

for i =1: length(mexFILE)

  mexcuda( mexcuda_opts{:}, sprintf('%s/%s.cu',mexPATH,mexFILE{i}), inc{1}, inc{2});

  system(sprintf('mv %s.mexa64 %s/compiled',mexFILE{i}, mexPATH));
end
end





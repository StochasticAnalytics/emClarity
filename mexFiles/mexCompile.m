function [] = mexCompile()

fprintf("\n\nCompile here\n\n");

mexPATH = '~/work/emClarity/mexFiles/';
% mexFILE = 'mexFFT';
mexFILE = {'mexCTF','mexFFT'};

mexcuda_opts = { ...
'-lcublas'          ...            % Link to cuBLAS
'-lmwlapack'        ...            % Link to LAPACK
'-lcufft'           ...            % Link to cuFFT
['NVCCFLAGS= -Wno-deprecated-gpu-targets']...wh
'-L/groups/grigorieff/home/himesb/thirdParty/cuda-10.0/lib64'   ...    % Location of CUDA libraries
'-L/groups/grigorieff/home/himesb/thirdParty/cuda-10.0/nvvm/lib64'};

% '-L/usr/local/cuda-9.1/lib64'   ...    % Location of CUDA libraries
% '-L/usr/local/cuda-9.1/nvvm/lib64'};
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/lib64'   ...    % Location of CUDA libraries
% '-L/groups/grigorieff/home/himesb/thirdParty/cuda-8.0/nvvm/lib64'};

for i =1: length(mexFILE)
  mexcuda(mexcuda_opts{:}, sprintf('%s%s.cu',mexPATH,mexFILE{i}));


  system(sprintf('mv %s.mexa64 %s/compiled',mexFILE{i}, mexPATH));
end
end

%nvcc -ptx --library-path /groups/grigorieff/home/himesb/thirdParty/cuda-9.2/nvvm/lib64



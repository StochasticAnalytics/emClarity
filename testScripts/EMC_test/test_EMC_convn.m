function tests = test_EMC_convn
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_convn;
testCase.TestData.evaluateOutput = @evaluateOutput;
testCase.TestData.fixtureRand = @help_getInputRand;
testCase.TestData.fixtureKernel = @EMC_gaussianKernel;
end


function [result, message] = evaluateOutput(IMAGE, KERNEL, OUTPUTCELL, EXTRA)
% Check that precision, method and size match the IMAGE
% If EXTRA (sigma), compute the standard kernel with EMC_gaussianKernel
% and check that output is equal to output with separable kernels. 

result = 'failed';

if length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    convolvedImg = OUTPUTCELL{1};
end

% precision and method
[actualMethod, actualPrecision] = help_getClass(convolvedImg);
[expectMethod, expectPrecision] = help_getClass(IMAGE);

if ~strcmp(actualMethod, expectMethod)
    message = sprintf('expected method=%s, got %s', expectMethod, actualMethod);
    return
elseif ~strcmp(actualPrecision, expectPrecision)
    message = sprintf('expected precision=%s, got %s', expectPrecision, actualPrecision);
    return
end

% size
if ~isequal(size(convolvedImg), size(IMAGE))
    message = sprintf('expected size=%s, got %s', size(IMAGE), size(convolvedImg));
    return
end

% separable
if EXTRA && isvector(KERNEL)
    kSize = zeros(1, ndims(IMAGE)) + length(KERNEL);
   	out = EMC_convn(IMAGE, ...
                    EMC_gaussianKernel(kSize, EXTRA, actualMethod, {'precision', actualPrecision}));
  	if any(abs(out - convolvedImg) > 1e-5, 'all')
      	message = 'separable kernel does not give the same results as standard kernel';
       	return
    end
end

result = 'passed';
message = '';

end


function test_default_2d(testCase)
% Test for 2d, cpu/gpu, single/double.
%
sigma = {1};
sizesKernel = {[5,5];[6,6];[1,5];[5,1];[1,6];[6,1]};

% cpu double
sizes = help_getRandomSizes(2, [10,300], '2d');
img = help_getBatch({'fixtureRand'}, {'cpu'}, {'double'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'cpu'}, {{'precision', 'double'}});
testCase.TestData.toTest = help_getBatch(img, kernel, {false}, sigma);

% cpu single
sizes = help_getRandomSizes(2, [10,300], '2d');
img = help_getBatch({'fixtureRand'}, {'cpu'}, {'single'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'cpu'}, {{'precision', 'single'}});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, kernel, {false}, sigma)];

% gpu double
sizes = help_getRandomSizes(2, [10,300], '2d');
img = help_getBatch({'fixtureRand'}, {'gpu'}, {'double'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'gpu'}, {{'precision', 'double'}});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, kernel, {false}, sigma)];

% cpu single
sizes = help_getRandomSizes(2, [10,300], '2d');
img = help_getBatch({'fixtureRand'}, {'gpu'}, {'single'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'gpu'}, {{'precision', 'single'}});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, kernel, {false}, sigma)];

EMC_runTest(testCase);

end


function test_default_3d(testCase)
% Test for 3d, cpu/gpu, single/double.
%
sigma = {1};
sizesKernel = {[5,5,5];[6,6,6];[1,5];[5,1];[1,6];[6,1]};

% cpu double
sizes = help_getRandomSizes(2, [10,300], '3d');
img = help_getBatch({'fixtureRand'}, {'cpu'}, {'double'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'cpu'}, {{'precision', 'double'}});
testCase.TestData.toTest = help_getBatch(img, kernel, {false}, sigma);

% cpu single
sizes = help_getRandomSizes(2, [10,300], '3d');
img = help_getBatch({'fixtureRand'}, {'cpu'}, {'single'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'cpu'}, {{'precision', 'single'}});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, kernel, {false}, sigma)];

% gpu double
sizes = help_getRandomSizes(2, [10,300], '3d');
img = help_getBatch({'fixtureRand'}, {'gpu'}, {'double'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'gpu'}, {{'precision', 'double'}});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, kernel, {false}, sigma)];

% cpu single
sizes = help_getRandomSizes(2, [10,300], '3d');
img = help_getBatch({'fixtureRand'}, {'gpu'}, {'single'}, sizes);
kernel = help_getBatch({'fixtureKernel'}, sizesKernel, sigma, {'gpu'}, {{'precision', 'single'}});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, kernel, {false}, sigma)];

EMC_runTest(testCase);
end


function test_assumption(testCase)

img = {ones(10,10); 12314; rand(1,10); rand(10,1); {}; {'wr', '124'}; ...
       nan; inf; ones(3,3,3,3); ones(0,10); ones(0,0); 'efqrg'; "wqgtg"};
kernel = {{}; {'wr', '124'}; 'efqrg'; "wqgtg"};
testCase.TestData.toTest = help_getBatch(img, kernel, {'EMC:IMAGE'}, {false});

% precision and method
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    {ones(10,10), ones(10,10, 'single'), 'EMC:IMAGE', false; ...
     ones(10,10, 'single'), ones(10,10), 'EMC:IMAGE', false; ...
     ones(10,10,'gpuArray'), ones(10,10, 'single','gpuArray'), 'EMC:IMAGE', false; ...
     ones(10,10, 'single','gpuArray'), ones(10,10,'gpuArray'), 'EMC:IMAGE', false; ...
     ones(10,10,'single','gpuArray'), ones(10,10, 'single'), 'EMC:IMAGE', false; ...
     ones(10,10,'gpuArray'), ones(10,10), 'EMC:IMAGE', false; ...
    }];
EMC_runTest(testCase);

end

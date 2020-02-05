function tests = test_EMC_helpFunctions

% Test for:
%   EMC_getClass
%   EMC_setMethod
%   EMC_setPrecision
%   EMC_isOnGpu
%   EMC_is3d
%   EMC_getOption
%

tests = functiontests(localfunctions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_EMC_getClass(testCase)
% [PRECISION, flgGPU, METHOD] = EMC_getClass(IMAGE)

testCase.TestData.functionToTest = @EMC_getClass;
testCase.TestData.evaluateOutput = @evaluate_getClass;
testCase.TestData.debug = 0;
testCase.TestData.testName = '';

testCase.TestData.toTest = { ...
    randn(100,100, 'single'),              false, {'single', 'cpu'}; ...
    randn(100,100),                        false, {'double', 'cpu'}; ...
    gpuArray(randn(100,100, 'single')),    false, {'single', 'gpu'}; ...
    randn(100,100, 'gpuArray'),            false, {'double', 'gpu'}; ...
    };

EMC_runTest(testCase);

end

function [result, message] = evaluate_getClass(~, OUTPUTCELL, EXTRA)
%
% Check the precision, method and flgGpu match the EXTRA.

result = 'passed';
message = '';

% For EMC_getClass, 3 outputs is expected.
if length(OUTPUTCELL) ~= 3
    result = 'failed';
    message = sprintf('expected output number = 3, got %d', length(OUTPUTCELL));
    return
else
    [precision, flgGpu, method] = OUTPUTCELL{:};
end

if ~(isscalar(flgGpu) && islogical(flgGpu))
    result = 'failed';
    message = sprintf('flgGpu should be a boolean, got %s, size: %s', class(flgGpu), mat2str(flgGpu));
    return
elseif ~(strcmp(method, 'cpu') || strcmp(method, 'gpu'))
    result = 'failed';
    message = sprintf('method should be cpu or gpu, got %s', method);
    return
elseif ~(strcmp(precision, 'single') || strcmp(precision, 'double'))
    result = 'failed';
    message = sprintf('precision should be single or double, got %s', precision);
    return
end

[expectedPrecision, expectedMethod] = EXTRA{:};

if ~strcmp(precision, expectedPrecision)
    result = 'failed';
    message = sprint('precision should be %s, got %s', expectedPrecision, precision);
elseif ~strcmp(method, expectedMethod)
    result = 'failed';
    message = sprint('method should be %s, got %s', expectedMethod, method);
elseif flgGpu && ~strcmp(method, 'gpu')
    result = 'failed';
    message = 'If method=gpu, flgGpu should be true, got false';
elseif ~flgGpu && ~strcmp(method, 'cpu')
    result = 'failed';
    message = 'If method=cpu, flgGpu should be false, got true';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_EMC_setMethod(testCase)
% NUM = EMC_setMethod(NUM, METHOD)

testCase.TestData.functionToTest = @EMC_setMethod;
testCase.TestData.evaluateOutput = @evaluate_setMethod;
testCase.TestData.debug = 0;
testCase.TestData.testName = '';

img = {randn(100,100, 'single'); ...
       randn(100,100); ...
       gpuArray(randn(100,100, 'single')); ...
       randn(100,100, 'gpuArray')};

method = {'cpu'; 'gpu'};
expectedError = {false};
extra = {false};

testCase.TestData.toTest = help_getBatch(img, method, expectedError, extra);
EMC_runTest(testCase);

end


function [result, message] = evaluate_setMethod(NUM, METHOD, OUTPUTCELL, ~)

result = 'passed';
message = '';

% For EMC_setMethod, 1 outputs is expected.
if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    outNum = OUTPUTCELL{:};
end

[actualMethod, actualPrecison] = help_getClass(outNum);
[~, expectedPrecison] = help_getClass(NUM);

if ~strcmp(actualMethod, METHOD)
    result = 'failed';
    message = sprintf('output method should be %s, got %s', METHOD, actualMethod);
elseif ~strcmp(expectedPrecison, actualPrecison)
    result = 'failed';
    message = 'precision was changed';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_setPrecision(testCase)
% NUM = EMC_setPrecision(NUM, PRECISION)

testCase.TestData.functionToTest = @EMC_setPrecision;
testCase.TestData.evaluateOutput = @evaluate_setPrecision;
testCase.TestData.debug = 0;
testCase.TestData.testName = '';

img = {randn(100,100, 'single'); ...
       randn(100,100); ...
       gpuArray(randn(100,100, 'single')); ...
       randn(100,100, 'gpuArray')};

precision = {'single'; 'double'};
expectedError = {false};
extra = {false};

testCase.TestData.toTest = help_getBatch(img, precision, expectedError, extra);
EMC_runTest(testCase);

end


function [result, message] = evaluate_setPrecision(NUM, PRECISION, OUTPUTCELL, ~)

result = 'passed';
message = '';

% For EMC_setPrecision, 1 outputs is expected.
if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    outNum = OUTPUTCELL{:};
end

[actualMethod, actualPrecison] = help_getClass(outNum);
[expectedMethod, ~] = help_getClass(NUM);

if ~strcmp(actualPrecison, PRECISION)
    result = 'failed';
    message = sprintf('output precision should be %s, got %s', PRECISION, actualPrecison);
elseif ~strcmp(expectedMethod, actualMethod)
    result = 'failed';
    message = 'method was changed';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_is3d(testCase)
% [is3d, SIZE, ndim] = EMC_is3d(SIZE)

testCase.TestData.functionToTest = @EMC_is3d;
testCase.TestData.evaluateOutput = @evaluate_is3d;
testCase.TestData.debug = 0;

testCase.TestData.toTest = { ...
    [10,10,10],	false, {true,   [10,10,10], 3}; ...
    [10,10,1],  false, {false,  [10,10],    2}; ...
    [1,10,10],  false, {true,   [1,10,10],  3}; ...
    [10,10],  	false, {false,  [10,10],    2}; ...
    
    % vectors
    [1, 10],    false, {false,  [1,10],     2}; ...
    [10, 1],    false, {false,  [10,1],     2}; ...
    
    % scalar are accepted
    [1, 1],     false, {false,  [1,1],      2}; ...
    
    
    [10,-10,10],	'EMC:SIZE',	false; ...
    [10,-10],       'EMC:SIZE',	false; ...
    [10,0,1],       'EMC:SIZE',	false; ...
    [1,10,1,1],     'EMC:SIZE',	false; ...
    {},             'EMC:SIZE',	false; ...
    ones(10,10),    'EMC:SIZE',	false; ...
    100,            'EMC:SIZE',	false; ...
    
    % only integers
    [10,10.1,10],	'EMC:SIZE',	false; ...
    [10.5, 1],      'EMC:SIZE',	false; ...
    
    [nan, 1],       'EMC:SIZE', false; ...
    [inf, 12],      'EMC:SIZE', false; ...
    nan,            'EMC:SIZE', false; ...
    inf,         	'EMC:SIZE', false; ...
    
    [10;10]         'EMC:SIZE', false; ...
    };

EMC_runTest(testCase);

end


function [result, message] = evaluate_is3d(~, OUTPUTCELL, EXTRA)
% [is3d, SIZE, ndim] = EMC_is3d(SIZE)

result = 'passed';
message = '';

% For EMC_is3d, 3 outputs are expected.
if length(OUTPUTCELL) ~= 3
    result = 'failed';
    message = sprintf('expected output number = 3, got %d', length(OUTPUTCELL));
    return
elseif ~isequal(OUTPUTCELL, EXTRA)
    result = 'failed';
    message = sprinft('outputs are not equal to the expected output in EXTRA');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_isOnGpu(testCase)
% ISONGPU = EMC_isOnGpu(NUM)

testCase.TestData.functionToTest = @EMC_isOnGpu;
testCase.TestData.evaluateOutput = @evaluate_isOnGpu;
testCase.TestData.debug = 0;

testCase.TestData.toTest = { ...
    randn(100,100, 'single'),           false, false; ...
    randn(100,100)                      false, false; ...
    gpuArray(randn(100,100, 'single')), false, true; ...
    randn(100,100, 'gpuArray'),         false, true; ...
    };

testCase.TestData.testName = '';
EMC_runTest(testCase);

end

function [result, message] = evaluate_isOnGpu(~, OUTPUTCELL, EXTRA)
% ISONGPU = EMC_isOnGpu(NUM)

result = 'passed';
message = '';

% For EMC_isOnGpu, 1 outputs is expected.
if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    isOnGpu = OUTPUTCELL{:};
end

if isOnGpu ~= EXTRA
    result = 'failed';
    message = sprinft('should be %d, got %d', isOnGpu, EXTRA);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_getOption(testCase)
% [OPTION] = EMC_getOption(OPTION, ONLY, FILTER)

testCase.TestData.functionToTest = @EMC_getOption;
testCase.TestData.evaluateOutput = @evaluate_getOption;
testCase.TestData.debug = 2;

inCell0 = {};
inCell1 = {'param1', 1};
inCell2 = {'param1', 1; 'param2', 2};
inCell3 = {'param1', 1; 'param2', 2; 'param3', 3};

inStru0 = struct;
inStru1 = cell2struct(inCell1(:, 2), inCell1(:, 1), 1);
inStru2 = cell2struct(inCell2(:, 2), inCell2(:, 1), 1);
inStru3 = cell2struct(inCell3(:, 2), inCell3(:, 1), 1);

testCase.TestData.toTest = {
    inCell0, {},                    false,  false,        	inStru0; ...
    inCell0, {'param1'},            true,   false,        	inStru0; ...
    inCell1, {'param2'},            false, 	'EMC:OPTION',   false; ...
    inCell1, {'param2'},            true, 	false,         	inStru0; ...
    inCell2, {'param1', 'param2'},  true, 	false,         	inStru2; ...
    inCell2, {'param1'},            true,  	false,         	inStru1; ...
    inCell2, {'param1', 'param2'},  true, 	false,         	inStru2; ...
    inCell3, {},                    false,  'EMC:OPTION',   false; ...
    inCell3, {},                    true,   false,        	inStru0; ...
    inCell3, {'param1', 'param3'},  true,   false,       	struct('param1', 1, 'param3', 3); ...
    
    inStru0, {},                    false,  false,        	inStru0; ...
    inStru0, {'param1'},            true,   false,       	inStru0; ...
    inStru1, {'param2'},            false, 	'EMC:OPTION',   false; ...
    inStru1, {'param2'},            true, 	false,         	inStru0; ...
    inStru2, {'param1', 'param2'},  true, 	false,        	inStru2; ...
    inStru2, {'param1'},            true,  	false,         	inStru1; ...
    inStru2, {'param1', 'param2'},  true, 	false,         	inStru2; ...
    inStru3, {},                    false,  'EMC:OPTION',   false; ...
    inStru3, {},                    true,   false,         	inStru0; ...
    inStru3, {'param1', 'param3'},  true,   false,         	struct('param1', 1, 'param3', 3); ...
    
    % assume error
    1, {},              false, 	'error', false; ...
    1:3, {},            false, 	'error', false; ...
    'kayak', {},       	false, 	'error', false; ...
    {'kayak', 1}, {},   false,	'error', false; ...
    {1, 1}, {},         false, 	'error', false; ...
    {{}, 1}, {},        false, 	'error', false; ...
    {{'kayak'}, 1}, {}, false,	'error', false; ...
    {1:3, 1}, {},       false, 	'error', false; ...
    {true, 1}, {},      false, 	'error', false; ...

    };

EMC_runTest(testCase);

end

function [result, message] = evaluate_getOption(~, ~, ~, OUTPUTCELL, EXTRA)
% [OPTION] = EMC_getOption(OPTION, ONLY, FILTER)

result = 'passed';
message = '';

% For EMC_getOption, 1 outputs is expected.
if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    output = OUTPUTCELL{:};
    if ~isstruct(output)
        result = 'failed';
        message = sprintf('output should be a structure, got %s', class(output));
        return
    end
end

if ~isequaln(output, EXTRA)
    result = 'failed';
    message = 'output structure is not equal to the expected structure';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_sharePrecision(testCase)
% SHARE = EMC_sharePrecision(NUM1, NUM2)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_sharePrecision;
testCase.TestData.evaluateOutput = @evaluate_sharePrecision;

testCase.TestData.toTest = {
    ones(10,10,'single'), ones(10,10,'single'), false, true; ...
    ones(10,10,'double'), ones(10,10,'double'), false, true; ...
    ones(10,10,'single'), ones(10,10,'double'), false, false; ...
    ones(10,10,'double'), ones(10,10,'single'), false, false; ...
    
    ones(10,10,'single','gpuArray'), ones(10,10,'single','gpuArray'), false, true; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'double','gpuArray'), false, true; ...
    ones(10,10,'single','gpuArray'), ones(10,10,'double','gpuArray'), false, false; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'single','gpuArray'), false, false; ...
    
    ones(10,10,'single','gpuArray'), ones(10,10,'single'), false, true; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'double'), false, true; ...
    ones(10,10,'single','gpuArray'), ones(10,10,'double'), false, false; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'single'), false, false; ...
    
    ones(10,10,'single'), ones(10,10,'single','gpuArray'), false, true; ...
    ones(10,10,'double'), ones(10,10,'double','gpuArray'), false, true; ...
    ones(10,10,'single'), ones(10,10,'double','gpuArray'), false, false; ...
    ones(10,10,'double'), ones(10,10,'single','gpuArray'), false, false; ...
    };

EMC_runTest(testCase);
end


function [result, message] = evaluate_sharePrecision(~, ~, OUTPUTCELL, EXTRA)
% SHARE = EMC_sharePrecision(NUM1, NUM2)

result = 'passed';
message = '';

if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    share = OUTPUTCELL{1};
end

if share ~= EXTRA
    result = 'failed';
    message = sprinft('should be %d, got %d', share, EXTRA);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_shareMethod(testCase)
% SHARE = EMC_shareMethod(NUM1, NUM2)

testCase.TestData.functionToTest = @EMC_shareMethod;
testCase.TestData.evaluateOutput = @evaluate_shareMethod;

testCase.TestData.toTest = {
    ones(10,10,'single'), ones(10,10,'single'), false, true; ...
    ones(10,10,'double'), ones(10,10,'double'), false, true; ...
    ones(10,10,'single'), ones(10,10,'double'), false, true; ...
    ones(10,10,'double'), ones(10,10,'single'), false, true; ...
    
    ones(10,10,'single','gpuArray'), ones(10,10,'single','gpuArray'), false, true; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'double','gpuArray'), false, true; ...
    ones(10,10,'single','gpuArray'), ones(10,10,'double','gpuArray'), false, true; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'single','gpuArray'), false, true; ...
    
    ones(10,10,'single','gpuArray'), ones(10,10,'single'), false, false; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'double'), false, false; ...
    ones(10,10,'single','gpuArray'), ones(10,10,'double'), false, false; ...
    ones(10,10,'double','gpuArray'), ones(10,10,'single'), false, false; ...
    
    ones(10,10,'single'), ones(10,10,'single','gpuArray'), false, false; ...
    ones(10,10,'double'), ones(10,10,'double','gpuArray'), false, false; ...
    ones(10,10,'single'), ones(10,10,'double','gpuArray'), false, false; ...
    ones(10,10,'double'), ones(10,10,'single','gpuArray'), false, false; ...
    
    {}, 1:10, false, true; ...
    };
EMC_runTest(testCase);
end


function [result, message] = evaluate_shareMethod(~, ~, OUTPUTCELL, EXTRA)
% SHARE = EMC_sharePrecision(NUM1, NUM2)

result = 'passed';
message = '';

if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    share = OUTPUTCELL{1};
end

if share ~= EXTRA
    result = 'failed';
    message = sprinft('should be %d, got %d', share, EXTRA);
end
end

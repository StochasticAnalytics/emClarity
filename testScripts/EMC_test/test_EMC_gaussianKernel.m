function tests = test_EMC_gaussianKernel
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_gaussianKernel;
testCase.TestData.evaluateOutput = @evaluateOutput;
end


function [result, message] = evaluateOutput(SIZE, ~, OPTION, OUTPUTCELL, ~)
% [KERNEL] = EMC_gaussianKernel(SIZE, SIGMA, OPTION)
% Check:
%   METHOD:     output has the expected method
%   PRECISION:  output has the expected precision
%   SIZE:       output has the desired size
%   The sum should be one.
%

result = 'failed';

if length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    kernel = OUTPUTCELL{1};
end

% Correct precision and method
% By default, should be single and cpu. Otherwise, respect OPTION parameters.
[actualMethod, actualPrecision] = help_getClass(kernel);

if iscell(OPTION) && ~isempty(OPTION) && any(strcmp(OPTION(:, 1), 'precision'))
    expectedPrecision = OPTION{strcmp(OPTION(:, 1), 'precision'), 2};
else
    expectedPrecision = 'single';  % default
end
if iscell(OPTION) && ~isempty(OPTION) && any(strcmp(OPTION(:, 1), 'method'))
    expectedMethod = OPTION{strcmp(OPTION(:, 1), 'method'), 2};
else
    expectedMethod = 'cpu';  % default
end

if ~strcmp(actualMethod, expectedMethod)
   	message = sprintf('expected output method = %s, got %s', actualMethod, expectedMethod);
    return
elseif ~strcmp(actualPrecision, expectedPrecision)
   	message = sprintf('expected output precision = %s, got %s', actualPrecision, expectedPrecision);
    return
end

% size
if ~isequal(SIZE, size(kernel))
    message = sprintf('expected size = %s, got %s', SIZE, size(kernel));
    return
end

% sum
if abs(sum(kernel(:)) - 1) > 1e-5
    message = sprintf('expected size = %s, got %s', SIZE, size(kernel));
    return
end

result = 'passed';
message = '';
end


function test_default(testCase)
% 2d
sizes = help_getRandomSizes(5, [3, 20], '2d');
sigmas = [help_getRandomSizes(2, [1,5], '2d'); help_getRandomSizes(2, [1,5], '1d')];
option = help_getBatchOption({'precision', {'single'; 'double'}; ...
                              'method', {'cpu'; 'gpu'}});

testCase.TestData.toTest = help_getBatch(sizes, sigmas, option, {false}, {false});

% 3d
sizes = help_getRandomSizes(2, [3, 20], '3d');
sigmas = [help_getRandomSizes(2, [1,5], '3d'); help_getRandomSizes(2, [1,5], '1d')];
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(sizes, sigmas, option, {false}, {false})];

EMC_runTest(testCase);
end


function test_assumption(testCase)
% sizes
sizes = {[5,5,5,5]; 5; 'kayak'; []; {}; nan; inf; [nan, nan]; [inf; inf]; [1,1,1]; [1,1]};
testCase.TestData.toTest = help_getBatch(sizes, {1}, {{}}, {'EMC:SIZE'}, {false});

% sigmas
sigmas = {[1,1,1]; 'kayak'; nan; inf; {}; -2; 0; [nan;1]};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, sigmas, {{}}, {'EMC:SIGMA'}, {false}); ...
                            {[10,10,10], [10,10], {}, 'EMC:SIGMA', false}];

% options
options = help_getBatchOption({'precision', {''; 'doubles'; 12; []; nan; inf}; ...
                               'method', {''; 'doubles'; 12; []; nan; inf}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {1}, options, {'error'}, {false})];
EMC_runTest(testCase);
end

function tests = test_EMC_gaussianKernel
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_gaussianKernel;
testCase.TestData.evaluateOutput = @evaluateOutput;
end


function [result, message] = evaluateOutput(SIZE, ~, METHOD, OPTION, OUTPUTCELL, ~)
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

if ~strcmpi(actualMethod, METHOD)
   	message = sprintf('expected output method = %s, got %s', actualMethod, METHOD);
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
method = {'cpu'; 'gpu'};
sigmas = [help_getRandomSizes(2, [1,5], '2d'); help_getRandomSizes(2, [1,5], '1d')];
option = help_getBatchOption({'precision', {'single'; 'double'}});

testCase.TestData.toTest = help_getBatch(sizes, sigmas, method, option, {false}, {false});

% 3d
sizes = help_getRandomSizes(2, [3, 20], '3d');
sigmas = [help_getRandomSizes(2, [1,5], '3d'); help_getRandomSizes(2, [1,5], '1d')];
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(sizes, sigmas, method, option, {false}, {false})];

EMC_runTest(testCase);
end


function test_assumption(testCase)
% sizes
sizes = {[5,5,5,5]; 5; 'kayak'; []; {}; nan; inf; [nan, nan]; [inf; inf]; [1,1,1]; [1,1]};
testCase.TestData.toTest = help_getBatch(sizes, {1}, {'cpu'}, {{}}, {'EMC:SIZE'}, {false});

% sigmas
sigmas = {[1,1,1]; 'kayak'; nan; inf; {}; -2; 0; [nan;1]};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, sigmas, {'cpu'}, {{}}, {'EMC:SIGMA'}, {false}); ...
                            {[10,10,10], [10,10], 'cpu', {}, 'EMC:SIGMA', false}];

% method
methods = {''; 'doubles'; 12; []; nan; inf};
testCase.TestData.toTest = help_getBatch({[10,10]}, {1}, methods, {{}}, {'EMC:METHOD'}, {false});

% options
options = help_getBatchOption({'precision', {''; 'doubles'; 12; []; nan; inf}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {1}, {'cpu'}, options, {'error'}, {false})];
EMC_runTest(testCase);
end

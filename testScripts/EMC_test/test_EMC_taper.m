function tests = test_EMC_taper

% Errors that can be raised by EMC_taper:
%

tests = functiontests(localfunctions);

end


function setupOnce(testCase)
% Setup
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_taper;
testCase.TestData.evaluateOutput = @evaluateOutput;

end


function [result, message] = evaluateOutput(~, SIZE, OPTION, OUTPUTCELL, ~)
% Output taper should be correct size and row numeric vector.
% The output precision and method depends on FIRST and LAST.

result = 'passed';
message = '';

% For EMC_taper, only one output is expected.
if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    OUTPUT = OUTPUTCELL{1};
end

% numeric vector
if ~isnumeric(OUTPUT) && ~isvector(OUTPUT)
    result = 'failed';
    message = sprintf('taper should be numeric vector, got class: %s, size: %s', ...
                      mat2str(class(OUTPUT)), mat2str(size(OUTPUT)));
    return
end

% correct length, column vector
sizeTaper = size(OUTPUT);
if sizeTaper(1) ~= 1 && SIZE ~= sizeTaper(2)
 	result = 'failed';
    message = sprintf('taper should be a row vector of %d elements, got size:%s', mat2str(size(OUTPUT)));
    return
end

% Correct precision and method
% By default, should be single and cpu. Otherwise, respect OPTION parameters.
[actualMethod, actualPrecision] = help_getClass(OUTPUT);

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
    result = 'failed';
   	message = sprintf('expected output method = %s, got %s', actualMethod, expectedMethod);
    return
elseif ~strcmp(actualPrecision, expectedPrecision)
    result = 'failed';
   	message = sprintf('expected output precision = %s, got %s', actualPrecision, expectedPrecision);
    return
end

% Last value should be 0 or OPTION.end
if iscell(OPTION) && ~isempty(OPTION) && any(strcmp(OPTION(:, 1), 'end'))
    endValue = OPTION{strcmp(OPTION(:, 1), 'end'), 2};
else
    endValue = 0;  % default
end
if endValue ~= OUTPUT(1, end)
    result = 'failed';
   	message = sprintf('last value should be %s, got %s', endValue, OUTPUT(1, end));
    return
end

% If first=true (include first value), first value should be 1 or OPTION.start
if iscell(OPTION) && ~isempty(OPTION) && any(strcmp(OPTION(:, 1), 'first'))
    if OPTION{strcmp(OPTION(:, 1), 'first'), 2}
        if iscell(OPTION) && ~isempty(OPTION) && any(strcmp(OPTION(:, 1), 'start'))
            startValue = OPTION{strcmp(OPTION(:, 1), 'start'), 2};
        else
            startValue = 1;  % default
        end
        if startValue ~= OUTPUT(1, 1)
            result = 'failed';
            message = sprintf('first value should be %d, got %d', startValue, OUTPUT(1, 1));
            return
        end
    end
end
    
end  % evaluateOutput


function test_default(testCase)

type_ = {'cosine'; 'linear'};
sizes = {100};
option = help_getBatchOption({'precision', {'single'; 'double'}; ...
                              'start', {1; -5; 2.5; single(1.5); gpuArray(2)}; ...
                              'end', {1; -5; 2.5; single(1.5); gpuArray(2)}; ...
                              'first', {true; false}; ...
                              'method', {'cpu'; 'gpu'}});
expectedError = {false};
extra = {false};

testCase.TestData.toTest = help_getBatch(type_, sizes, option, expectedError, extra);
EMC_runTest(testCase);

end


function test_assumption(testCase)
%% type
type_ = {[]; ''; {}; 'kayak'; 2; 1:12};
sizes = {100};
option = {{}};

expectedError = {'EMC_taper:TYPE'};
extra = {false};

testCase.TestData.toTest = help_getBatch(type_, sizes, option, expectedError, extra);
testCase.TestData.toTest = [testCase.TestData.toTest; {{'cosine'}, 1, {}, 'EMC_taper:SIZE', false}];

%% optinal parameters
type_ = {'cosine'};
sizes = {10};
option = help_getBatchOption({'precision', {''; 'kayak'; 1; false; []; {12}}; ...
                              'method', {''; 'kayak'; 1; false; {}; {12}};
                              'start', {1:3; 'kayak'; {}}; ...
                              'end', {1:3; ''; {}}; ...
                              'first', {2; '1'; nan}});
option = option(~cellfun(@isempty, option));

expectedError = {'error'};  % generic error
extra = {false};

testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(type_, sizes, option, expectedError, extra)];
EMC_runTest(testCase);

end

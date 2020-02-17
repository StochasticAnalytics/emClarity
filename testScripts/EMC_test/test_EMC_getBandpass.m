function tests = test_EMC_getBandpass
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_getBandpass;
testCase.TestData.evaluateOutput = @evaluateOutput;
end


function [result, message] = evaluateOutput(SIZE, ~, ~, ~, METHOD, OPTION, OUTPUTCELL, EXTRA)
%
% Check:
%   METHOD:     output is following the METHOD instruction.
%   PRECISION:  output has the desired precision
%   SIZE:       output has the correct size.
%
%   If EXTRA, load and check for equality with fixture.
%

result = 'failed';

if length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    bandpass = OUTPUTCELL{1};
end

% precision and method
[actualMethod, actualPrecision] = help_getClass(bandpass);
if help_isOptionDefined(OPTION, 'precision')
    expectedPrecision = help_getOptionParam(OPTION, 'precision');
else
    expectedPrecision = 'single';  % default
end
if ~strcmp(actualMethod, METHOD)
    message = sprintf('expected method=%s, got %s', METHOD, actualMethod);
    return
elseif ~strcmp(actualPrecision, expectedPrecision)
    message = sprintf('expected precision=%s, got %s', expectedPrecision, actualPrecision);
    return
end
    
% size
expectedSize = size(bandpass);
if SIZE ~= expectedSize
 	message = sprintf('expected size=%d, got %d', expectedSize, SIZE);
   	return
end
    
% compare with fixture
if EXTRA
    if any(bandpass - EXTRA) > 1e-7
       	message = 'vectors are not equal to corresponding fixture';
       	return
    end
end

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default(testCase)
% Test for 2d, cpu/gpu and every option.
%

sizes = [help_getRandomSizes(1, [50, 1000], '2d'); help_getRandomSizes(1, [50, 250], '3d')];
pixels = {1; 2.5};
highpass = {5; nan; 0};
lowpass = {5; nan; 0; 'nyquist'};
method = {'cpu'; 'gpu'};

option = help_getBatchOption({'origin', {-1; 0; 1; 2}; ...
                              'precision', {'single' ; 'double'}; ...
                              'half', {true; false}});

testCase.TestData.toTest = help_getBatch(sizes, pixels, highpass, lowpass, ...
                                         method, option, {false}, {false});
EMC_runTest(testCase);

end


function test_roll(testCase)

sizes = help_getRandomSizes(1, [50, 1000], '2d');
pixels = {2.5};
highpass = {30};
lowpass = {10};
method = {'cpu'};

option = help_getBatchOption({'origin', {-1; 0}; ...
                              'half', {false}; ...
                              'highpassRoll', {0.05; 0.1; 1; 'extended'}; ...
                              'lowpassRoll', {0.05; 0.1; 1; 'extended'}; ...
                              'highpassThresh', {0; 0.1; 0.99}; ...
                              'lowpassThresh', {0; 0.1; 0.99}});

testCase.TestData.toTest = help_getBatch(sizes, pixels, highpass, lowpass, ...
                                         method, option, {false}, {false});
EMC_runTest(testCase);
                          
end


function test_assumptions(testCase)
% make sure error are raised

% sizes
sizes = {'kayak'; [1,10]; [10,1]; [10,10,10,10]; []; 2; {}; [0,10]; [10.5;10]; [-10,10]; [10;10]};
testCase.TestData.toTest = help_getBatch(sizes, {1}, {1}, {1}, {'cpu'}, {{}}, {'EMC:SIZE'}, {false});

% pixel
pixels = {'kayak'; -1; {}; [1,2]; []; nan; inf; 0};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch([100,100], pixels, {1}, {1}, {'cpu'}, {{}}, {'EMC:PIXEL'}, {false})];

% highpass
highpass = {'kayak'; -1; {}; [1,2]; []; inf};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch([100,100], {1}, highpass, {2}, {'cpu'}, {{}}, {'EMC:HIGHPASS'}, {false})];

% lowpass
lowpass = {'kayak'; -1; {}; [1,2]; []; inf; 0.5};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch([100,100], {1}, {2}, lowpass, {'cpu'}, {{}}, {'EMC:LOWPASS'}, {false})];

% method
methods = {'kayak'; 1; {}; [1,2]; []; nan; inf};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch([100,100], {1}, {2}, {2}, methods, {{}}, {'EMC:METHOD'}, {false})];

% option
options = help_getBatchOption({'origin', {-12; 3; 'kayak'; nan; inf; {'231'}}; ...
                               'half', {1;0;-1;inf;nan;[];''}; ...
                               'precision', {'int'; 'float'; 12; []; nan; inf}});
options = [options; help_getBatchOption({'highpassRoll', {''; 'kayak'; [1,2]; {}}; ...
                                         'lowpassRoll', {''; 'kayak'; [1,2]; {}}; ...
                                         'highpassThresh', {-0.1; 1; nan; inf; ''; [1,2]}; ...
                                         'lowpassThresh', {-0.1; 1; nan; inf; ''; [1,2]}})];
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch([100,100], {1}, {2}, {2.5}, {'cpu'}, options, {'error'}, {false})];

EMC_runTest(testCase);

end


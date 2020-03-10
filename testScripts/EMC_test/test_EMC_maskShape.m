function tests = test_EMC_maskShape
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_maskShape;
testCase.TestData.evaluateOutput = @evaluateOutput;
end


function [result, message] = evaluateOutput(~, SIZE, ~, METHOD, OPTION, OUTPUTCELL, EXTRA)
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
    mask = OUTPUTCELL{1};
end

% precision and method
[actualMethod, actualPrecision] = help_getClass(mask);
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
expectedSize = size(mask);
if SIZE ~= expectedSize
 	message = sprintf('expected size=%d, got %d', expectedSize, SIZE);
   	return
end
    
% compare with fixture
if EXTRA
    if any(mask - EXTRA) > 1e-7
       	message = 'vectors are not equal to corresponding fixture';
       	return
    end
end

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default_2d(testCase)
% Test for 2d, cpu/gpu and every option.
%

shapes = {'rectangle'; 'sphere'; 'cylinder'};
sizes = help_getRandomSizes(1, [50, 1000], '2d');
radius = [help_getRandomSizes(1, [10, 400], '2d'); [30.5,32.1]];
method = {'cpu'; 'gpu'};

kernel = [1.5e-06 0.00013 0.0044 0.054 0.24 0.4 0.24 0.054 0.0044 0.00013 1.5e-06];
option = help_getBatchOption({'origin', {1; 2}; ...
                              'precision', {'single' ; 'double'}; ...
                              'sym', {1;2;6}; ...
                              'kernel', {true;false;0;0.2;kernel}; ...
                              'shift', [help_getRandomSizes(1, [-100,100], '2d'); [30.5,-32.1]]});

testCase.TestData.toTest = help_getBatch(shapes, sizes, radius, method, option, {false}, {false});
EMC_runTest(testCase);

end


function test_default_3d(testCase)
% Test for 3d, cpu/gpu and every option.
%

shapes = {'rectangle'; 'sphere'; 'cylinder'};
sizes = help_getRandomSizes(1, [50, 200], '3d');
radius = [help_getRandomSizes(1, [10, 100], '3d'); [30.5,32.1,20.4]];
method = {'cpu'; 'gpu'};

kernel = [1.5e-06 0.00013 0.0044 0.054 0.24 0.4 0.24 0.054 0.0044 0.00013 1.5e-06];
option = help_getBatchOption({'origin', {1; 2}; ...
                              'precision', {'single' ; 'double'}; ...
                              'sym', {1;2;6}; ...
                              'kernel', {true;false;0;0.2;kernel}; ...
                              'shift', [help_getRandomSizes(1, [-100,100], '3d'); [-30.5,32.1,-20.4]]});

testCase.TestData.toTest = help_getBatch(shapes, sizes, radius, method, option, {false}, {false});
EMC_runTest(testCase);

end


function test_default_outOfRange(testCase)
% Test weird radius and shifts: is there an out of range error?
% test only 3d (2d is the same).

shapes = {'rectangle'; 'sphere'; 'cylinder'};
sizes = help_getRandomSizes(1, [150, 200], '3d');
radius = help_getRandomSizes(1, [150, 300], '3d');
method = {'cpu'; 'gpu'};

option = help_getBatchOption({'origin', {1;2}; ...
                              'sym', {2;6}; ...
                              'kernel', {false}; ...
                              'shift', help_getRandomSizes(3, [-200,200], '3d')});

testCase.TestData.toTest = help_getBatch(shapes, sizes, radius, method, option, {false}, {false});
EMC_runTest(testCase);

end


function test_assumptions(testCase)
% make sure error are raised

% shape
shapes = {''; 'spheres'; 1; {}; {'2'}; nan; inf};
testCase.TestData.toTest = help_getBatch(shapes, {[20,20]}, {[20,20]}, {'cpu'}, {{}}, {'EMC:SHAPE'}, {false});

% sizes
sizes = {'sphere'; [1,10]; [10,1]; [10,10,10,10]; []; 2; {}; [0,10]; [10.5;10]; ...
         [-10,10]; [nan, 10]; [inf,10]; [10;10]};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, sizes, {[10,10]}, {'cpu'}, {{}}, {'EMC:SIZE'}, {false})];

% radius
radius = {'sphere'; [10,10,10,10]; []; 2; {}; [0,10]; ...
          [10.5;10]; [-10,10]; [nan, 10]; [inf,10]; [10;10]};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, radius, {'cpu'}, {{}}, {'EMC:RADIUS'}, {false}); ...
    {'rectangle', [10,10,10], [10,10], 'cpu', {}, 'EMC:RADIUS', false}];

% method
methods = {'kayak'; 1; {}; [1,2]; []; nan; inf};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, methods, {{}}, {'EMC:METHOD'}, {false})];

% options
options = help_getBatchOption({'origin', {-1;3;'kayak'; nan; inf; {}}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, {'cpu'}, options, {'error'}, {false})];

options = help_getBatchOption({'sym', {1.2; 0; -1; 'kayak'; nan; inf; {}}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, {'cpu'}, options, {'error'}, {false})];

options = help_getBatchOption({'kernel', {1.1; -0.3; 12; nan; inf; ''; []}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, {'cpu'}, options, {'error'}, {false})];

options = help_getBatchOption({'shift', {1.1; -0.3; 12; nan; inf; ''; [1,1,1]; [10;10]}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, {'cpu'}, options, {'error'}, {false})];

options = help_getBatchOption({'precision', {'int'; 'float'; 12; []; nan; inf}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, {'cpu'}, options, {'error'}, {false})];

testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'sphere'}, {[10,10]}, {[10,10]}, {'cpu'}, options, {'error'}, {false}); ...
    {'rectangle', [10,10,10], [10,10,10], 'cpu', {'shift', [10,10]}, 'EMC:shift', false; ...
     'rectangle', [10,20], [10,10], 'cpu', {}, 'EMC:kernel', false}];

EMC_runTest(testCase);

end

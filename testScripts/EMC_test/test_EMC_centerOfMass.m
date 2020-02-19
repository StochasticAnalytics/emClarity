function tests = test_EMC_centerOfMass
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_centerOfMass;
testCase.TestData.evaluateOutput = @evaluateOutput;
testCase.TestData.fixtureRand = @help_getInputRand;
end


function [result, message] = evaluateOutput(IMAGE, ~, OUTPUTCELL, EXTRA)
% COM = EMC_centerOfMass(IMAGE, ORIGIN)

result = 'failed';

if length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    COM = OUTPUTCELL{1};
end

% precision and method
[actualMethod, actualPrecision] = help_getClass(COM);
[expectMethod, expectPrecision] = help_getClass(IMAGE);

if ~strcmp(actualMethod, expectMethod)
    message = sprintf('expected method=%s, got %s', expectMethod, actualMethod);
    return
elseif ~strcmp(actualPrecision, expectPrecision)
    message = sprintf('expected precision=%s, got %s', expectPrecision, actualPrecision);
    return
end

% size
if size(COM) ~= ndims(IMAGE)
 	message = sprintf('expected size=%d, got %d', ndims(IMAGE), size(COM));
   	return
end

% compare with fixture
if EXTRA
    if any(abs(COM - EXTRA) > 1e-7, 'all')
       	message = 'com is not equal to corresponding fixture';
       	return
    end
end

result = 'passed';
message = '';

end


function test_default(testCase)
% COM = EMC_centerOfMass(IMAGE, ORIGIN)
%

sizes = help_getRandomSizes(2, [500, 1000], '2d');
precisions = {'single'; 'double'};
methods = {'cpu'; 'gpu'};
img = help_getBatch({'fixtureRand'}, methods, precisions, sizes);

origins = {0;1;2};

testCase.TestData.toTest = help_getBatch(img, origins, {false}, {false});

sizes = help_getRandomSizes(2, [100,300], '3d');
img = help_getBatch({'fixtureRand'}, methods, precisions, sizes);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, origins, {false}, {false})];

EMC_runTest(testCase);

end


function test_assumptions(testCase)
% make sure error are raised

% image
imgNan = ones(10,10);
imgNan(4,4) = nan;
imgInf = imgNan;
imgInf(4,4) = inf;
img = {'kayak'; ones(5,5,5,5); []; 2; {}; imgNan; imgInf};
testCase.TestData.toTest = help_getBatch(img, {1}, {'error'}, {false});

% origin
origin = {'kayak'; ones(1,10); ones(10,1); ones(5,5,5,5); []; -1; 0.3; nan; inf; {}};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({ones(10,10)}, origin, {'EMC:origin'}, {false})];
EMC_runTest(testCase);

end


function test_fixture(testCase)
testCase.TestData.toTest = { ...
    ones(10,10), 1, false, [-0.5,-0.5]; ...
    ones(10,10), 2, false, [0.5,0.5]; ...
    ones(10,10), 0, false, [0,0]};

EMC_runTest(testCase);

end
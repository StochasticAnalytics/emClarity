function tests = perf_EMC_resize
tests = functiontests(localfunctions);
end


function setupOnce(testCase)

% Setup
testCase.TestData.functionToTest = @EMC_resize;
testCase.TestData.fixtureOnes = @help_getInputOnes;
end


function test_without_taper(testCase)

%% 2d
inSizes = {[100,100]; [1000,1000]; [4000,4000]; [8000,8000]};
precision = {'single'; 'double'};

limits = {[0,0,0,0]; [30,30,30,30]; [-30,-30,-30,-30]; [30,-30,-30,30]};
option = {{'taper', false}};

img = help_getBatch({'fixtureOnes'}, {'cpu'}, precision, inSizes);
testCase.TestData.toTest = help_getBatch(img, limits, option, {'cpu'});

img = help_getBatch({'fixtureOnes'}, {'gpu'}, precision, inSizes);
testCase.TestData.toTest = [testCase.TestData.toTest; help_getBatch(img, limits, option, {'gpu'})];
testCase.TestData.testName = '2d';
EMC_runPerf(testCase);

%% 3d
inSizes = {[50,50,50]; [200,200,200]};
precision = {'single'; 'double'};

limits = {[0,0,0,0,0,0]; [15,15,15,15,15,15]; [-15,-15,-15,-15,-15,-15]; [15,-15,15,-15,15,-15]};
option = {{'taper', false}};

img = help_getBatch({'fixtureOnes'}, {'cpu'}, precision, inSizes);
testCase.TestData.toTest = help_getBatch(img, limits, option, {'cpu'});

img = help_getBatch({'fixtureOnes'}, {'gpu'}, precision, inSizes);
testCase.TestData.toTest = [testCase.TestData.toTest; help_getBatch(img, limits, option, {'gpu'})];
testCase.TestData.testName = '3d';

EMC_runPerf(testCase);
end



function test_with_taper(testCase)

%% 2d
inSizes = {[100,100]; [1000,1000]; [4000,4000]; [8000,8000]};
precision = {'single'; 'double'};

limits = {[0,0,0,0]; [30,30,30,30]; [30,30,0,0]; [-30,-30,-30,-30]; [30,-30,-30,30]};
option = {{'taper', true}};

img = help_getBatch({'fixtureOnes'}, {'cpu'}, precision, inSizes);
testCase.TestData.toTest = help_getBatch(img, limits, option, {'cpu'});

img = help_getBatch({'fixtureOnes'}, {'gpu'}, precision, inSizes);
testCase.TestData.toTest = [testCase.TestData.toTest; help_getBatch(img, limits, option, {'gpu'})];
testCase.TestData.testName = '2d';
EMC_runPerf(testCase);

%% 3d
inSizes = {[50,50,50]; [200,200,200]};
precision = {'single'; 'double'};

limits = {[0,0,0,0,0,0]; [15,15,15,15,15,15]; [15,15,0,0,0,0]; [-15,-15,-15,-15,-15,-15]; [15,-15,15,-15,15,-15]};
option = {{'taper', true}};

img = help_getBatch({'fixtureOnes'}, {'cpu'}, precision, inSizes);
testCase.TestData.toTest = help_getBatch(img, limits, option, {'cpu'});

img = help_getBatch({'fixtureOnes'}, {'gpu'}, precision, inSizes);
testCase.TestData.toTest = [testCase.TestData.toTest; help_getBatch(img, limits, option, {'gpu'})];
testCase.TestData.testName = '3d';

EMC_runPerf(testCase);
end


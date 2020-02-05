function tests = test_EMC_resize
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_resize;
testCase.TestData.evaluateOutput = @evaluateOutput;
testCase.TestData.fixtureOnes = @help_getInputOnes;
end


function [result, message] = evaluateOutput(IMG, LIMITS, OPTION, OUTPUTCELL, ~)
% This function is checking for output size, precision and method.

result = 'passed';
message = '';

% For EMC_resize, only one output is expected.
if length(OUTPUTCELL) ~= 1
    result = 'failed';
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    OUTPUT = OUTPUTCELL{1};
end

outSize = size(IMG) + sum(reshape(LIMITS, 2, []));
if numel(outSize) == 3 && outSize(3) == 1
    outSize = outSize(1:2);
end

[actualMethod, actualPrecision] = help_getClass(OUTPUT);
[expectedMethod, expectedPrecision] = help_getClass(IMG);
if iscell(OPTION) && ~isempty(OPTION) && any(strcmp(OPTION(:, 1), 'precision'))
    expectedPrecision = OPTION{strcmp(OPTION(:, 1), 'precision'), 2};
end

if ~strcmp(actualMethod, expectedMethod)
    result = 'failed';
    message = sprintf('img method was changed from %s to %s', expectedMethod, actualMethod);
elseif ~strcmp(actualPrecision, expectedPrecision)
    result = 'failed';
    message = sprintf('expected precision: %s, got %s', expectedPrecision, actualPrecision);
elseif outSize ~= size(OUTPUT)
    result = 'failed';
    message = sprintf('output size should be %s, got %s', mat2str(outSize), mat2str(size(OUTPUT)));
end

end


function [result, message] = evaluateOutputFixture(IMG, LIMITS, OPTION, OUTPUTCELL, EXTRA)
% This function is checking for output size, precision and method.
% Additionnaly, checks that the output is equal to fixture.

[result, message] = evaluateOutput(IMG, LIMITS, OPTION, OUTPUTCELL, nan);
if strcmp(result, 'failed')
    return
end

fixture = load(EXTRA);
if any(abs(OUTPUTCELL{1} - fixture.out)) > 1e-07
    result = 'failed';
    message = 'not equal to fixture';
end

end


function test_default(testCase)
%
% 1) output size matches the input size + limits. Many combinaison are tried, as well as
%    specific cases (crop/pad only, no crop/pad, both, one side only, one axis only)
% 2) For each run, the precision of the input image should stay unchanged by default.
%    Single and double precision are tested.
% 3) Make sure method stays unchanged. Both cpu and gpu are tested.
%

%% 2d
imgSizes = help_getRandomSizes(1, [80, 1000], '2d');
img = help_getBatch({'fixtureOnes'}, {'cpu'; 'gpu'}, {'single'; 'double'}, imgSizes);

limits = help_getRandomLimits(1, [-30,30], {'both', 'pad', 'crop'}, ...
                              {'all', 'x', 'y'}, {'both', 'left', 'right'}, '2d');

option = {{}; {'origin', -1}; {'origin', 1}; {'origin', 2}};

testCase.TestData.toTest = help_getBatch(img, limits, option, {false}, {false});
testCase.TestData.testName = '2d';
EMC_runTest(testCase);

%% 3d
imgSizes = help_getRandomSizes(1, [80, 250], '3d');
img = help_getBatch({'fixtureOnes'}, {'cpu'; 'gpu'}, {'single'; 'double'}, imgSizes);

limits = help_getRandomLimits(1, [-30,30], {'both', 'pad', 'crop'}, ...
                              {'all', 'x', 'y', 'z', 'xy', 'xz', 'yz'}, {'both', 'left', 'right'}, '3d');

option = {{}; {'origin', -1}; {'origin', 1}; {'origin', 2}};

testCase.TestData.toTest = help_getBatch(img, limits, option, {false}, {false});
testCase.TestData.testName = '3d';
EMC_runTest(testCase);

end % test_default


function test_precision(testCase)
% 1) output size matches the input size + limits.
% 2) precision is respected and changed if whished.
% 3) method stays unchanged.
% 4) 3d is correctly handled.
% 5) reciprocal space is correctly handled.
% 6) This test make sure the precision of the taper doesn't influence the output precision.
% 7) padding values are correctly cast to the desired precision, as well as 'uniform' mode.
% 8) errors are correctly raised: only 'single' and 'double' precision are allowed.

%% 2d
imgSizes = help_getRandomSizes(1, [80, 1000], '2d');
img = help_getBatch({'fixtureOnes'}, {'cpu'; 'gpu'}, {'single'; 'double'}, imgSizes);

limits = help_getRandomLimits(2, [-30,30], {'both'}, {'all',}, {'both'}, '2d');

option = help_getBatchOption({'origin', {-1;1}; ...
                              'taper', {single([0.7,0.5,0.3,0.1]); single(2.5)}; ...
                              'value', {2.5; single(2.4); 'uniform'; 'mean'}});

testCase.TestData.toTest = help_getBatch(img, limits, option, {false}, {false});
testCase.TestData.testName = 'batch2d';
EMC_runTest(testCase);

%% 3d
imgSizes = help_getRandomSizes(1, [80, 200], '3d');
img = help_getBatch({'fixtureOnes'}, {'cpu'; 'gpu'}, {'single'; 'double'}, imgSizes);

limits = help_getRandomLimits(2, [-30,30], {'both'}, {'all',}, {'both'}, '3d');


testCase.TestData.toTest = help_getBatch(img, limits, option, {false}, {false});
testCase.TestData.testName = 'batch3d';
EMC_runTest(testCase);

%% assume errors

option = help_getBatchOption({'precision', {1; 1.4; 'kayak'; []; [1,2,3]; {}; nan; inf}});
option = option(~cellfun(@isempty, option));
testCase.TestData.toTest = help_getBatch(img, limits, option, {'EMC:precision'}, {false});
testCase.TestData.testName = 'assumeError';
EMC_runTest(testCase);

end % test_precision


function test_taper_and_value(testCase)

%% 2d
imgSizes = help_getRandomSizes(1, [120, 1000], '2d');
img = help_getBatch({'fixtureOnes'}, {'cpu'; 'gpu'}, {'single'; 'double'}, imgSizes);

limits = help_getRandomLimits(2, [-30,30], {'both'}, {'all',}, {'both'}, '2d');

option = help_getBatchOption({'origin', {-1;1}; ...
                              'taper', {single([0.7,0.5,0.3,0.1]); single(2.5); 5; ...
                                        gpuArray(single([0.7,0.5,0.3,0.1])); ...
                                        gpuArray(single(2.5)); {'cosine', 10}; {'linear', 30}}; ...
                              'value', {single(2.4); gpuArray(single(2.4)); gpuArray(12); 'uniform'; 'mean'};
                              'force_taper', {true}});

testCase.TestData.toTest = help_getBatch(img, limits, option, {false}, {false});
testCase.TestData.testName = 'batch2d';
EMC_runTest(testCase);

%% 3d
imgSizes = help_getRandomSizes(1, [120, 200], '3d');
img = help_getBatch({'fixtureOnes'}, {'cpu'; 'gpu'}, {'single'; 'double'}, imgSizes);

limits = help_getRandomLimits(2, [-30,30], {'both'}, {'all',}, {'both'}, '3d');

option = help_getBatchOption({'origin', {-1;1}; ...
                              'taper', {single([0.7,0.5,0.3,0.1]); single(2.5); 5; ...
                                        gpuArray(single([0.7,0.5,0.3,0.1])); ...
                                        gpuArray(single(2.5)); {'cosine', 10}; {'linear', 30}}; ...
                              'value', {single(2.4); gpuArray(single(2.4)); gpuArray(12); 'uniform'; 'mean'};
                              'force_taper', {true}});

testCase.TestData.toTest = help_getBatch(img, limits, option, {false}, {false});
testCase.TestData.testName = 'batch3d';
EMC_runTest(testCase);

%% assume errors
% Origin and 2d/3d does not matter as an error should be raised within the checkIN.

% EMC:taper
option = help_getBatchOption({'taper', {[0.9; 0.7; 0.5; 0.3; 0.1]; [0.9, 0.7, 0.5; 0.3, 0.1, 0]; ...
                                        'kayak'; []; {}; ""; struct('cosine', 10)}});
option = option(~cellfun(@isempty, option));
testCase.TestData.toTest = help_getBatch(img, limits, option, {'EMC:taper'}, {false});

% EMC:taper - SIZE
option = help_getBatchOption({'taper', {{'cosine', 1}; {'cosine', -1}; {'linear', 1}; {'linear', -1}; ...
                                        {'cosine', nan}; {'cosine', inf}}});
option = option(~cellfun(@isempty, option));
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, limits, option, {'EMC:taper'}, {false})];

% EMC:taper - TYPE
option = help_getBatchOption({'taper', {{'kayak', 10}; {[], 10}; {nan, 10}; {inf, 10}}});
option = option(~cellfun(@isempty, option));
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, limits, option, {'EMC:taper'}, {false})];
                        
% EMC:IMAGE
img1 = {'fixtureOnes', 'cpu', 'single', [100,100]};
img2 = {'fixtureOnes', 'cpu', 'single', [100,101]};

testCase.TestData.toTest = [testCase.TestData.toTest; { ...
    img1, [10,10,10,10], {'taper', {'linear', 101}}, 'EMC:taper', false; ...
    img2, [10,10,10,10], {'taper', {'linear', 102}}, 'EMC:taper', false; ...
    img2, [10,10,10,10], {'taper', {'linear', 51}; 'origin', -1}, 'EMC:taper', false; ...
    }];

% EMC:value
option = help_getBatchOption({'value', {{'kayak'}; []; ""; 'kayak'; [1,2]}});
option = option(~cellfun(@isempty, option));
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(img, limits, option, {'EMC:value'}, {false})];

testCase.TestData.testName = 'assumeError';
EMC_runTest(testCase);

end  % test_taper


function test_assumptions(testCase)
img = {'fixtureOnes', 'cpu', 'single', [20,20]};

testCase.TestData.toTest = { ...
    % origin should be -1, 1 or 2.
    img, [5,5,5,5], {'origin', 0},                      'EMC:origin',    false; ...
    img, [5,5,5,5], {'origin', 3},                      'EMC:origin',    false; ...
    img, [5,5,5,5], {'origin', 1.2},                    'EMC:origin',    false; ...
    img, [5,5,5,5], {'origin', 'kayak'},                'EMC:origin',    false; ...
    img, [5,5,5,5], {'origin', nan},                    'EMC:origin',    false; ...
    img, [5,5,5,5], {'origin', inf},                    'EMC:origin',    false; ...
    
    % force taper should be a bool
    img, [5,5,5,5], {'force_taper', 3},                 'EMC:force_taper',   false; ...
    img, [5,5,5,5], {'force_taper', 1.2},               'EMC:force_taper',   false; ...
    img, [5,5,5,5], {'force_taper', 'kayak'},           'EMC:force_taper',   false; ...
    img, [5,5,5,5], {'force_taper', [true, false]'},    'EMC:force_taper',   false; ...
    img, [5,5,5,5], {'force_taper', nan},               'EMC:force_taper',   false; ...
    img, [5,5,5,5], {'force_taper', inf},               'EMC:force_taper',   false; ...
    
    % Img should be 2d or 3d numerical array
    ones(1,10),         [5,5,5,5], {},                 	'EMC:IMAGE',    false; ...
    10,                 [5,5,5,5], {},                 	'EMC:IMAGE',    false; ...
    'kayak',            [5,5,5,5], {},                 	'EMC:IMAGE',    false; ...
    {},                 [5,5,5,5], {},                  'EMC:IMAGE',   	false; ...
    1:-1,               [5,5,5,5], {},                  'EMC:IMAGE',   	false; ...
    nan,                [5,5,5,5], {},                  'EMC:IMAGE',   	false; ...
    inf,                [5,5,5,5], {},                  'EMC:IMAGE',  	false; ...
    
    [],                 [5,5,5,5], {},                  'EMC:SIZE',   	false; ...
    ones(10,10,10,10),  [5,5,5,5], {},                  'EMC:SIZE',  	false; ...
    
    % limits should be a numerical vector of size=ndim*2
    ones(20,20,20), [5,5,5,5], {},                      'EMC:LIMITS', false; ...
    img, [5,5;5,5;5,5], {},                             'EMC:LIMITS', false; ...
    img, [5,5,5,5,5],   {},                             'EMC:LIMITS', false; ...
    img, [5,5],         {},                             'EMC:LIMITS', false; ...
    img, [5,5;5,5],     {},                             'EMC:LIMITS', false; ...
    img, 5,             {},                             'EMC:LIMITS', false; ...
    img, 'kayak',       {},                             'EMC:LIMITS', false; ...
    img, [],            {},                             'EMC:LIMITS', false; ...
    img, {},            {},                             'EMC:LIMITS', false; ...
    img, nan,            {},                          	'EMC:LIMITS', false; ...
    img, inf,            {},                          	'EMC:LIMITS', false; ...
    img, [nan,10,10,10],            {},                 'EMC:LIMITS', false; ...
    img, [inf,10,10,10],            {},                 'EMC:LIMITS', false; ...
};

EMC_runTest(testCase);

end


function test_fixture(testCase)

testCase.TestData.evaluateOutput = @evaluateOutputFixture;

img2d = {'fixtureOnes', 'cpu', 'single', [128,128]};
img3d = {'fixtureOnes', 'cpu', 'single', [128,128,128]};
fixt = 'fixtures/EMC_resize_fixture_';

testCase.TestData.toTest = { ...
	img2d, [0,0,0,0],           {'taper', false},                                   false, [fixt, '001']; ...
    img2d, [0,0,0,0],           {'taper', false; 'origin', -1},                     false, [fixt, '001']; ...
    img3d, [0,0,0,0,0,0],       {'taper', false},                                   false, [fixt, '002']; ...
    img3d, [0,0,0,0,0,0],       {'taper', false; 'origin', -1},                     false, [fixt, '002']; ...
         
    img2d, [0,0,0,0],           {'force_taper', true},                              false, [fixt, '003']; ...
    img2d, [0,0,0,0],           {'force_taper', true; 'origin', -1},                false, [fixt, '004']; ...
    img3d, [0,0,0,0,0,0],       {'force_taper', true},                              false, [fixt, '005']; ...
    img3d, [0,0,0,0,0,0],       {'force_taper', true; 'origin', -1},                false, [fixt, '006']; ...
         
    img2d, [0,0,0,0],           {'force_taper', true; 'value', 5},                  false, [fixt, '007']; ...
    img2d, [0,0,0,0],           {'force_taper', true; 'origin', -1; 'value', 5},    false, [fixt, '008']; ...
    img3d, [0,0,0,0,0,0],       {'force_taper', true; 'value', 5},                  false, [fixt, '009']; ...
    img3d, [0,0,0,0,0,0],       {'force_taper', true; 'origin', -1; 'value', 5}, 	false, [fixt, '010']; ...
         
    img2d, [10,11,11,10],       {'taper', false},                                   false, [fixt, '011']; ...
    img2d, [10,11,11,10],       {'taper', false; 'origin', -1},                     false, [fixt, '012']; ...
    img3d, [10,11,11,10,5,6],	{'taper', false},                                   false, [fixt, '013']; ...
    img3d, [10,11,11,10,5,6],   {'taper', false; 'origin', -1},                     false, [fixt, '014']; ...
         
    img2d, [10,11,11,10],       {'taper', false; 'value', 5},                       false, [fixt, '015']; ...
    img2d, [10,11,11,10],       {'taper', false; 'origin', -1; 'value', 5},         false, [fixt, '016']; ...
    img3d, [10,11,11,10,5,6],   {'taper', false; 'value', 5},                       false, [fixt, '017']; ...
    img3d, [10,11,11,10,5,6],   {'taper', false; 'origin', -1; 'value', 5},         false, [fixt, '018']; ...
    };

EMC_runTest(testCase);

end

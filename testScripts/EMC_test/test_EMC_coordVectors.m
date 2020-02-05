function tests = test_EMC_coordVectors
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_coordVectors;
testCase.TestData.evaluateOutput = @evaluateOutput;
end


function [result, message] = evaluateOutput(SIZE, METHOD, OPTION, OUTPUTCELL, EXTRA)
%
% Check:
%   METHOD:     outputs are following the METHOD instruction.
%   PRECISION:  outputs have the desired precision
%   SIZE:       outputs have the correct size.
%               If unit dimension (size=1), check that output is NaN.
%
%   If EXTRA, load and check for equality with fixture.
%

result = 'failed';
ndim = numel(SIZE);
if ndim == 3 && SIZE(3) == 1
    ndim = 2;
end

if length(OUTPUTCELL) ~= 3
    message = sprintf('expected output number = 3, got %d', length(OUTPUTCELL));
    return
elseif ndim == 2 && (~isscalar(OUTPUTCELL{3}) && ~isnan(OUTPUTCELL{3}))
    message = sprintf('for ndim=%d, vZ should be NaN', ndim);
    return
end

for iVec = 1:ndim
    vector = OUTPUTCELL{iVec};
    
    % nan
    if SIZE(iVec) == 1
        if ~isscalar(vector) && ~isnan(vector)
            message = sprintf('vector %d should be nan', iVec);
            return
        end
        continue
    end
    
    % precision and method
    [actualMethod, actualPrecision] = help_getClass(vector);
    
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
    expectedSize = SIZE(iVec);
    
    % half = true
    if iVec == 1 && help_isOptionDefined(OPTION, 'half')
        half = help_getOptionParam(OPTION, 'half');
        if isscalar(half) && islogical(half) && half  % half = true
            if help_isOptionDefined(OPTION, 'origin')
                origin = help_getOptionParam(OPTION, 'origin');
            else
                origin = 1;
            end
            
            if mod(SIZE(1), 2) || origin ~= 0  % odd or not real center
                expectedSize = floor(SIZE(1)/2) + 1;
            else  % even with origin=0
                expectedSize = SIZE(1)/2;
            end
        end
    end

    if numel(vector) ~= expectedSize
        message = sprintf('expected size=%d, got %d', expectedSize, numel(vector));
        return
    end
    
    % row vector
    if ~isrow(vector)
        message = sprintf('expected row vector, got size:%s', mat2str(size(vector)));
        return
    end
end

% compare with fixture
if iscell(EXTRA) && ~isempty(EXTRA)
    for iVec = 1:length(EXTRA)
        if any(abs(OUTPUTCELL{iVec} - EXTRA{iVec}) > 1e-7)
            message = 'vectors are not equal to corresponding fixture';
            return
        end
    end
end

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default(testCase)
% Test for 3d/2d/1d, cpu/gpu and every option except 'half'.
%

sizes = {[1, 30, 1]; [30, 1, 1]; [30, 30, 1]; [1, 30]; [30, 1]; [1, 1, 1]; [1, 1, 30]};
sizes = [sizes; help_getRandomSizes(1, [200, 400], '2d')];
sizes = [sizes; help_getRandomSizes(1, [200, 400], '3d')];

method = {'cpu'; 'gpu'};
option = help_getBatchOption({'origin', {0; 1; 2; -1}; ...
                              'normalize', {true; false}; ...
                              'isotrope', {true; false}; ...
                              'precision', {'single' ; 'double'}; ...
                              'half', {true; false}});

testCase.TestData.toTest = help_getBatch(sizes, method, option, {false}, {false});
EMC_runTest(testCase);

end


function test_shift(testCase)
% Separate test because shifts are not compatible with origin=-1 and half=true.
%
testCase.TestData.debug = 0;

%% Batch
for iDim = {'3d', '2d'}
    sizes = help_getRandomSizes(5, [2, 100], iDim{:});
    method = {'cpu'; 'gpu'};
    option = help_getBatchOption({'origin', {0; 1; 2}; ...
                                  'normalize', {true; false}; ...
                                  'isotrope', {true; false}; ...
                                  'precision', {'single' ; 'double'}; ...
                                  'shift', help_getRandomSizes(3, [-200, 200], iDim{:})});

    testCase.TestData.toTest = help_getBatch(sizes, method, option, {false}, {false});
    testCase.TestData.testName = iDim{:};
    EMC_runTest(testCase);
end

%% Assume errors and specific cases.
testCase.TestData.toTest = { ...
    [50,51,52], 'cpu', {'origin', -1; 'shift', [1,1,1]},                'EMC:shift', false; ...
    [51,52,53], 'cpu', {'half', true; 'shift', [1,1,1]},                'EMC:shift', false; ...
    [50,51,52], 'cpu', {'origin', -1; 'shift', [1,1,1]; 'half', true}, 	'EMC:shift', false; ...
    [50,51,52], 'cpu', {'shift', [nan,1,1]},                            'EMC:shift', false; ...
    [50,51,52], 'cpu', {'shift', [inf,1,1]},                            'EMC:shift', false; ...
    [50,51,52], 'cpu', {'shift', 'kayak'},                              'EMC:shift', false; ...
    [50,51],    'cpu', {'shift', [1,1,1]},                              'EMC:shift', false; ...
    [50,51],    'cpu', {'shift', [1,1,0]},                              'EMC:shift', false; ...
    [50,51],    'cpu', {'shift', [1;1]},                                'EMC:shift', false; ...
    };

end


function test_assumeError(testCase)

testCase.TestData.debug = 2;

% size
sizes = {[]; 1; [0,10]; [10,0]; [nan, 10]; [inf, 10]; '23'; ones(10,10); {}; nan; inf; [10;10]};
testCase.TestData.toTest = help_getBatch(sizes, {'cpu'}, {{}}, 'EMC:SIZE', {false});

% method
method = {''; {}; ""; "kayak"; [1,2]; 1; nan; inf};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, method, {{}}, 'EMC:METHOD', {false})];

% option
option = help_getBatchOption({'origin', {{}; 3; -2; 1.1; [1,2]; true; 'kayak'; nan; inf}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {'cpu'}, option, 'EMC:origin', {false})];

option = help_getBatchOption({'normalize', {{}; 3; 1.1; [true, false]; 1; 0; 'kayak'; nan; inf}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {'cpu'}, option, 'EMC:normalize', {false})];
                        
option = help_getBatchOption({'isotrope', {{}; 3; -2; 1.1; [1,2]; 1; 0; 'kayak'; nan; inf}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {'cpu'}, option, 'EMC:isotrope', {false})];

option = help_getBatchOption({'precision', {'cpu'; ''; {}; ""; "kayak"; [1,2]; 1; nan; inf}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {'cpu'}, option, 'EMC:precision', {false})];
                        
option = help_getBatchOption({'half', {{}; 3; 1.1; [true, false]; 1; 0; 'kayak'; nan; inf}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {'cpu'}, option, 'EMC:half', {false})];

EMC_runTest(testCase);

end


function test_fixture(testCase)

% Fixture

% size=10
fixt001 = [-5 -4 -3 -2 -1 0 1 2 3 4];                       % size=10, origin=1
fixt002 = [-4 -3 -2 -1 0 1 2 3 4 5];                        % size=10, origin=2
fixt003 = [-4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5];   % size=10, origin=0
fixt004 = [0 1 2 3 4 -5 -4 -3 -2 -1];                       % size=10, origin=-1

% size=11
fixt005 = [-5 -4 -3 -2 -1 0 1 2 3 4 5];                     % size=11, origin=1
fixt006 = [-5 -4 -3 -2 -1 0 1 2 3 4 5];                     % size=11, origin=2
fixt007 = [-5 -4 -3 -2 -1 0 1 2 3 4 5];                     % size=11, origin=0
fixt008 = [0 1 2 3 4 5 -5 -4 -3 -2 -1];                     % size=11, origin=-1

% size=10 with shifts
fixt009 = fixt001 - 2.5;                                    % size=10, origin=1
fixt010 = fixt002 + 4.5;                                    % size=10, origin=2
fixt011 = fixt003 - 0.5;                                    % size=10, origin=0

% size=11 with shifts
fixt012 = fixt005 + 2.5;                                    % size=11, origin=1
fixt013 = fixt006 - 11.5;                                   % size=11, origin=2
fixt014 = fixt007 - 0.5;                                    % size=11, origin=0

testCase.TestData.toTest = { ...
    [10,11], 'cpu', {},                                     false, {fixt001, fixt005}; ...
    [10,11], 'cpu', {'origin', 1},                          false, {fixt001, fixt005}; ...
    [10,11], 'cpu', {'origin', 2},                          false, {fixt002, fixt006}; ...
    [10,11], 'cpu', {'origin', 0},                          false, {fixt003, fixt007}; ...
    [10,11], 'cpu', {'origin',-1},                          false, {fixt004, fixt008}; ...
    
    [10,11], 'cpu', {'origin', 1; 'shift', [2.5, -2.5]},    false, {fixt009, fixt012}; ...
    [10,11], 'cpu', {'origin', 2; 'shift', [-4.5, 11.5]},	false, {fixt010, fixt013}; ...
    [10,11], 'cpu', {'origin', 0; 'shift', [0.5, 0.5]},     false, {fixt011, fixt014}; ...
};

EMC_runTest(testCase);

end

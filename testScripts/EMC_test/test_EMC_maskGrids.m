function tests = test_EMC_maskGrids
tests = functiontests(localfunctions);
end


function setupOnce(testCase)

% Setup
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_maskGrids;
testCase.TestData.evaluateOutput = @evaluateOutput;

end


function [result, message] = evaluateOutput(TYPE, SIZE, METHOD, OPTION, OUTPUTCELL, EXTRA)
%
% Check method, precision, size for each grid.
% If that that none of the size axis is 1 (execpt z), because EMC_maskGrids should raise an error for this.
% If EXTRA, load and check for equality with fixture.
%

result = 'failed';
ndim = numel(SIZE);
if ndim == 3 && SIZE(3) == 1
    ndim = 2;
end

if length(OUTPUTCELL) ~= 6
    message = sprintf('expected output number = 6, got %d', length(OUTPUTCELL));
    return
end

% Should have raised an error.
if any(SIZE(1:2) <= 1)
    message = sprintf('SIZE should not have empty or 1 pixel in X and Y, got %s', mat2str(SIZE));
    return
end

% The vectors are NOT checked and assumed to be correct.
% The vectors are tested by test_EMC_coordVectors.m
if ndim == 2
    expectedSize = [size(OUTPUTCELL{4}), size(OUTPUTCELL{5})];
else
    expectedSize = [size(OUTPUTCELL{4}), size(OUTPUTCELL{5}), size(OUTPUTCELL{6})];
end

for iGrid = 1:ndim
    grid = OUTPUTCELL{iGrid};
    
    % If radial, gY and gZ should be NaN.
    if strcmpi(TYPE, 'radial') && iGrid > 1
        if ~isscalar(grid) && ~isnan(grid)
            message = sprintf('grid %d should be nan', grid);
            return
        end
        continue
    end
    
    % precision and method
    [actualMethod, actualPrecision] = help_getClass(grid);
    
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
    if size(grid) ~= expectedSize
        message = sprintf('expected size=%d, got %d', expectedSize, size(grid));
        return
    end
end

% Compare for equality with fixture.
if isstring(EXTRA) || ischar(EXTRA)
    fixture = load(EXTRA);
    
    if strcmpi(TYPE, 'radial')
        if any(abs(OUTPUTCELL{1} - fixture.gX) > 1e-7)
            message = 'grids are not equal to corresponding fixture';
            return
        end
    else
        if any(abs(OUTPUTCELL{1} - fixture.gX) > 1e-7)
            message = 'gX is not equal to corresponding fixture';
            return
        end
        if any(abs(OUTPUTCELL{2} - fixture.gY) > 1e-7)
            message = 'gY is not equal to corresponding fixture';
            return
        end
        if any(abs(OUTPUTCELL{3} - fixture.gZ) > 1e-7)
            message = 'gZ is not equal to corresponding fixture';
            return
        end
    end
end

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default(testCase)

sizes = help_getRandomSizes(2, [1000, 4000], '2d');
sizes = [sizes; help_getRandomSizes(2, [50, 200], '3d')];

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
%% 2d
sizes = help_getRandomSizes(5, [2, 1000], '2d');
method = {'cpu'; 'gpu'};
option = help_getBatchOption({'origin', {0; 1; 2}; ...
                              'normalize', {true; false}; ...
                              'isotrope', {true; false}; ...
                              'precision', {'single' ; 'double'}; ...
                              'shift', help_getRandomSizes(3, [-200, 200], '2d')});

testCase.TestData.toTest = help_getBatch(sizes, method, option, {false}, {false});
testCase.TestData.testName = '2d';
EMC_runTest(testCase);


%% 3d
sizes = help_getRandomSizes(5, [2, 300], '3d');
method = {'cpu'; 'gpu'};
option = help_getBatchOption({'origin', {0; 1; 2}; ...
                              'normalize', {true; false}; ...
                              'isotrope', {true; false}; ...
                              'precision', {'single' ; 'double'}; ...
                              'shift', help_getRandomSizes(3, [-200, 200], '3d')});

testCase.TestData.toTest = help_getBatch(sizes, method, option, {false}, {false});
testCase.TestData.testName = '3d';
EMC_runTest(testCase);


end


function test_fixture(testCase)

testCase.TestData.toTest = { ...
    'radial', [30,51,30], 'cpu', {}, false,     './fixtures/EMC_maskGrids_fixture_001'; ...
    'radial', [30,51],    'cpu', {}, false,     './fixtures/EMC_maskGrids_fixture_002'; ...

    'cartesian', [30,51,30], 'cpu', {}, false,  './fixtures/EMC_maskGrids_fixture_003'; ...
    'cartesian', [30,51], 'cpu', {}, false,     './fixtures/EMC_maskGrids_fixture_004'; ...

    'spherical', [30,51,30], 'cpu', {}, false,  './fixtures/EMC_maskGrids_fixture_005'; ...
    'spherical', [30,51,30], 'cpu', {}, false,  './fixtures/EMC_maskGrids_fixture_006'; ...

    'cylinder', [30,51,30], 'cpu', {}, false,   './fixtures/EMC_maskGrids_fixture_007'; ...
    'cylinder', [30,51], 'cpu', {}, false,      './fixtures/EMC_maskGrids_fixture_008'; ...
    };

EMC_runTest(testCase);

end

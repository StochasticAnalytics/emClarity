function tests = test_EMC_limits
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_limits;
testCase.TestData.evaluateOutput = @evaluateOutput;
end


function [result, message] = evaluateOutput(CURRENT, DESIRED, OPTION, OUTPUTCELL, EXTRA)
% Check:
%   check limits respect desired size and preserve origin.
%   equal to fixture.

result = 'failed';
message = '';

% one output
if length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    actualLimits = reshape(OUTPUTCELL{1},2,[]);
end

% desired size is expected
outputSize = sum(actualLimits) + CURRENT;
if any(outputSize ~= DESIRED, 'all')
	message = sprintf('desired size is not respected: expected:%s, output:%s', ...
                      mat2str(DESIRED), mat2str(outputSize));
  	return
end  

% use another algorithm, more intuitive, to compute the limits and compare
expectedLimits = zeros(2, numel(CURRENT));
for iDim = 1:numel(CURRENT)
    current = CURRENT(iDim);
    desired = DESIRED(iDim);
    
    % start with
    if help_isOptionDefined(OPTION, 'origin')
        if abs(help_getOptionParam(OPTION, 'origin')) == 1  % 1 or -1
            if mod(current,2); left = 1; else; left = 0; end
        else
            if mod(current,2); left = 0; else; left = 1; end
        end
    else
      	if mod(current,2); left = 1; else; left = 0; end
    end
    
    difference = desired - current;
    % pad
    if difference > 0
        for i = 1:difference
            if left
                expectedLimits(1, iDim) = expectedLimits(1, iDim) + 1;
            else
                expectedLimits(2, iDim) = expectedLimits(2, iDim) + 1;
            end
            if left == 1; left = 0; else; left = 1; end  % switch side
        end
            
    % crop
    elseif difference < 0
        for i = 1:abs(difference)
            if left == 1; left = 0; else; left = 1; end  % switch side
            if left
                expectedLimits(1, iDim) = expectedLimits(1, iDim) - 1;
            else
                expectedLimits(2, iDim) = expectedLimits(2, iDim) - 1;
            end
            
        end
    end
end

% add shifts
if help_isOptionDefined(OPTION, 'shift')
	expectedLimits(1, :) = expectedLimits(1, :) + help_getOptionParam(OPTION, 'shift');
    expectedLimits(2, :) = expectedLimits(2, :) - help_getOptionParam(OPTION, 'shift');
end

if any(expectedLimits ~= actualLimits, 'all')
    message = sprintf('limits (%s) is not equal to expected limits (%s)', ...
                      mat2str(actualLimits), mat2str(expectedLimits));
    return
end

%% fixture
if ~islogical(EXTRA)
    if any(abs(limits - EXTRA) > 1e-7)
        message = 'limit is not equal to fixture';
        return
    end
end

result = 'passed';
end


function test_default(testCase)
% 2d
current = help_getRandomSizes(5, [1,500], '2d');
desired = help_getRandomSizes(5, [1,500], '2d');

option = help_getBatchOption({'origin', {1; -1; 2}});

testCase.TestData.toTest = help_getBatch(current, desired, option, {false}, {false});

% 3d
current = help_getRandomSizes(5, [1,500], '3d');
desired = help_getRandomSizes(5, [1,500], '3d');

option = help_getBatchOption({'origin', {1; -1; 2}});

testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(current, desired, option, {false}, {false})];
EMC_runTest(testCase);

end


function test_shift(testCase)
% origin=-1 is not allowed with shifts.
%

% 2d
current = help_getRandomSizes(3, [1,500], '2d');
desired = help_getRandomSizes(3, [1,500], '2d');

option = help_getBatchOption({'origin', {1; 2}; ...
                              'shift', help_getRandomSizes(3, [-100,100], '2d')});

testCase.TestData.toTest = help_getBatch(current, desired, option, {false}, {false});

% 3d
current = help_getRandomSizes(3, [1,500], '3d');
desired = help_getRandomSizes(3, [1,500], '3d');

option = help_getBatchOption({'origin', {1; 2}; ...
                              'shift', help_getRandomSizes(3, [-100,100], '3d')});

testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch(current, desired, option, {false}, {false})];
EMC_runTest(testCase);

end


function test_assumptions(testCase)

% current
current = help_getBatch({[nan, 10]; [10,nan]; inf; nan; [inf, 10]; ...
                         [10,inf]; [10.3,10]; [-9,10]; [0,20]; [10;10]});
testCase.TestData.toTest = help_getBatch(current, {[10,10]}, {{}}, {'EMC:LIMITS'}, {false});

% desired
desired = help_getBatch({[nan, 10]; [10,nan]; inf; nan; [inf, 10]; ...
                         [10,inf]; [10.3,10]; [-9,10]; [0,20]; [10;10]});
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, desired, {{}}, {'EMC:LIMITS'}, {false})];
                        
% options
option = help_getBatchOption({'origin', {3; 1.1; 0; [true, false]; -2; 'kayak'; nan; inf}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {[10,10]}, option, {'EMC:origin'}, {false})];

option = help_getBatchOption({'shift', {3; 1.1; -2; 'kayak'; nan; inf; [10;10]; [10,10,10]; [nan,10]; [inf,10]}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10]}, {[10,10]}, option, {'EMC:shift'}, {false})];

option = help_getBatchOption({'shift', {3; 1.1; -2; 'kayak'; nan; inf; [10,10]; ...
                                        [nan,10,10]; [inf,10,10]; [10;10;10]}});
option = option(~cellfun(@isempty, option), 1);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({[10,10,10]}, {[10,10,10]}, option, {'EMC:shift'}, {false}); ...
                            {[10,10,10], [10,10,10], {'shift', [10,10,10]; 'origin',-1}, 'EMC:shift', false}];

EMC_runTest(testCase);
end

function EMC_runPerf(testCase)
%
% Run perfomance tests (via timeit, gputimeit or your own function(s)).
%
% Inputs: testCase.TestData should contains the following field:
%   testName (string):          Name of the test function.
%
%   functionToTest (handle):    Function to test.
%                               Signature: [ouput1, ...] = example(input1, ...)
%                               Functions without outputs and/or inputs are also accepted.
%
%   toTest (cell):            	Should correspond to functionToTest inputs (same order), plus 1 arguments:
%                               - method (str|handle):
%                                   If str:     'cpu' or 'gpu'.
%                                   If handle:  function handle that will be called to perform the perfomance
%                                               measurement. This is usually useful if one wants to compare
%                                               gpu and cpu result.
%                                               Signature: [executionTime[int|float]] = funcHandle(input1, ...)
%
%                               dimension:      cell(nTests, (nInputs + 1))
%                               format:         {input1A, input1B, input1C, ..., method;
%                                                input2A, input2B, input2C, ..., method;
%                                                ...
%                                                inputNA, inputNB, inputNC, ..., method}
%
%   Optional field:
%     -> fixture_* (handle):    Same as EMC_runTest.
%
% Notes:
%   -
%
% Examples:
%   -
%
% Other m-files required:
%
% Created: 18Jan2020
% Last revision: 18Jan2020
% See also EMC_runTest
%

functionName = func2str(testCase.TestData.toTest);
logFileName = [pwd, '/logPerf/%s_%s', functionName, testCase.TestData.testName];
save(testCase.TestData.toTest, logFileName);

toTest = testCase.TestData.toTest;
executionTimes = cell(length(toTest), 1);

for iTest = 1:length(toTest)  % for every test
    
    % Extract fixture.
    nbArgs = length(testCase.TestData.toTest(iTest, :)) - 1;  % method
    cleanArg = cell(1, nbArgs);
    for iArg = 1:nbArgs  % for every input arguments
        currentArg = testCase.TestData.toTest{iTest, iArg};
        
        if (iscell(currentArg) && ~isempty(currentArg) && ...
            (isstring(currentArg{1}) || ischar(currentArg{1})) && ...
            contains(currentArg{1}, 'fixture_'))

            cleanArg(iArg) = testCase.TestData.toTest(iTest, iArg);
            if numel(currentArg) == 1
             	testCase.TestData.toTest{iTest, iArg} = testCase.TestData.(currentArg{1})();
            else
                testCase.TestData.toTest{iTest, iArg} = testCase.TestData.(currentArg{1})(currentArg{2:end});
            end
        end
    end
    
    % Run performation test.
    if isstring(toTest{iTest, end})
        if strcmpi(toTest{iTest, end}, 'cpu')
            executionTimes(iTest, 1) = timeit(@()testCase.TestData.functionToTest(toTest{iTest, 1:end-1}));
        elseif strcmpi(toTest{iTest, end}, 'gpu')
            executionTimes(iTest, 1) = gputimeit(@()testCase.TestData.functionToTest(toTest{iTest, 1:end-1}));
        else
            executionTimes(iTest, 1) = sprintf("method should be 'cpu', 'cpu' or function handle, got %s", ...
                                               toTest{iTest, end});
        end
    elseif isa(toTest{iTest, end}, 'function_handle')
        executionTimes(iTest, 1) = timeit(@()toTest{iTest, end}(toTest{iTest, 1:end-1}));
    else
        executionTimes(iTest, 1) = sprintf("method should be 'cpu', 'cpu' or function handle, got %s", ...
                                           toTest{iTest, end});
    end

    % Clean fixture.
    if any(~cellfun(@isempty, cleanArg))
        for iArg = 1:length(cleanArg)
            if ~isempty(cleanArg{iArg})
                testCase.TestData.toTest{iTest, iArg} = cleanArg{iArg};
            end
        end
    end
end

% Every tests were ran. Update the log with the results.
save([testCase.TestData.toTest, executionTimes], logFileName);

end  % EMC_runTest

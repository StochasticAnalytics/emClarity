function EMC_runTest(testCase)
%
% EMC functions testing framework.
%
% Objective:
%   1) Fast and robust way to implement tests for EMC functions.
%   2) Log: save test inputs, outputs and results into a mat file.
%   3) MATLAB debug: pause when a test fails.
%
% Input:
% testCase.TestData should contains the following field:
%
%   functionToTest (handle):    Function to test.
%                               Signature: [ouput1, ...] = example(input1, ...)
%
%   toTest (cell):              Should correspond to functionToTest inputs (same order), plus 2 arguments:
%                               - expectedError (false|str): If false: the test should not raise any error.
%                                                            If str: error identifier that should be raised.
%                                                                    If 'error', don't check for a specific
%                                                                    error id (generic error).
%
%                               - extraSpace (any): This variable is not used directly by EMC_runTest,
%                                                   but will be send to the evaluateOutput function.
%                                                   This allow more freedom when designing tests.
%
%                               dimensions: cell(nTests, (nInputs + 2))
%                               format: {input1A, input1B, input1C, ..., expectedError, extraSpace;
%                                        input2A, input2B, input2C, ..., expectedError, extraSpace;
%                                        ...
%                                        inputNA, inputNB, inputNC, ..., expectedError, extraSpace}
%
%   evaluateOutput (handle):    Function to check the validity of the output.
%                               It is only ran if the functionToTest did NOT raise an error.
%                               Signature: [result, message] = example(input1, ..., outputToCheck, extraSpace)
%                                   result (str):         whether or not the output pass the test;
%                                                         'passed', 'failed' or 'warning'.
%                                   message (str):    	  For logging purposes.
%                                   outputToCheck (cell): 1xn cell array containing the output to check.
%                                                         If the functionToTest does not return any output,
%                                                         the outputCell will be empty.
%                                   extraSpace (any):     See toTest field for more details.
%
%   (optional)
%   fixture* (handle):          If an input (input1A, etc.) is a 1xN cell with its first element starting by
%                               'fixture', EMC_runTest will look for the corresponding function handle in
%                               testCase.TestData, use the rest of the element(s) as input to run the
%                               function handle and use this output as input for the test. This mechanism
%                               allows to create or load fixture on the fly. Functions with no inputs are
%                               also accepted. See example below.
%
%   (optional)
%   debug (int):                0: no debug.
%                               1: if test result = 'fail', pause the execution.
%                               2: if test result = 'fail' or 'warning', pause the execution.
%
%   (optional)
%   testName (string):       	Name of the test function. Add this suffix to the logfilename.
%                               Use this prefix to have multiple calls within one test function.
%
% Notes:
%   - Logging:  Each test script (test_*.m) generates a directory called: ./logTest/log_<functionToTest>.
%               Within this directory, a log file is saved every time this function is called. The name of
%               this logfile is <testfunction><testNames>.m or <testfunction>.m if testName is not defined.
%
%   - Stack:    This function cannot be called from the cmd line. It needs a parent workspace.
%
%   - Add the possibility to test functions without outputs and/or inputs.
%
% Example:
%   % One test with EMC_resize:
%   testCase.TestData.functionToTest = @EMC_resize;
%
%   % Create inputs with one fixtures.
%   testCase.TestData.fixtureImg = @(Size, Precision) ones(Size, Precision);
%   img = {'fixtureImg', 'single', [128,128]};  % the array will be created just before running the test.
%   limits = [0,0,0,0];
%   option = {'taper', false};
%   expectedError = false;
%   extraSpace = nan;
%   testCase.TestData.toTest = {img, limits, option, expectedError, extraSpace};
%
%   % You should have created an evaluation function, with the following signature:
%   % [result, message] = exampleFunction(img, limits, option, outcell, extraSpace).
%   testCase.TestData.evaluateOutput = @exampleFunction;
%
%   % Run the test.
%   EMC_runTest(testCase);
%

% Created:	19Jan2020
% Version:  v.1.0.  new logging format (TF, 22Jan2020).
%           v.1.1.  debug is now optional (TF, 2Feb2020).
%

% Set debug pause.
if isfield(testCase.TestData, 'debug') && testCase.TestData.debug > 0
    dbstop if caught error EMC_runTest:debug
end

%% Save inputs as log.
stack = dbstack;
if length(stack) == 1  % EMC_runTest is called from the cmd line.
    error('This function should be called from a script, not from the command line.')
end
nameParentFunction = stack(2).name;
nameFunctionToTest = ['log_', func2str(testCase.TestData.functionToTest)];

folderName = [pwd, '/logTest/', nameFunctionToTest];
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

nbTests = length(testCase.TestData.toTest(:, 1));

if isfield(testCase.TestData, 'testName')
    logFileName = [folderName, sprintf('/%s_%s', nameParentFunction, testCase.TestData.testName)];
    fprintf('\t - %s:%s: %d tests\n', nameParentFunction, testCase.TestData.testName, nbTests)
else
    logFileName = [folderName, sprintf('/%s', nameParentFunction)];
    fprintf('\t - %s: %d tests\n', nameParentFunction, nbTests)
end

save(logFileName, 'testCase');

%% Start test
results = cell(nbTests, 2);  % result, message
for iTest = 1:nbTests
    
    % Extract fixture.
    nbArgs = length(testCase.TestData.toTest(iTest, :)) - 2;  % minus expectedError and extraSpace
    cleanArg = cell(1, nbArgs);
    for iArg = 1:nbArgs  % for every input arguments
        currentArg = testCase.TestData.toTest{iTest, iArg};
        
        if (iscell(currentArg) && ~isempty(currentArg) && ...
            (isstring(currentArg{1}) || ischar(currentArg{1})) && ...
            strncmpi(currentArg{1}, 'fixture', 7))

            cleanArg(iArg) = testCase.TestData.toTest(iTest, iArg);
            if numel(currentArg) == 1
             	testCase.TestData.toTest{iTest, iArg} = testCase.TestData.(currentArg{1})();
            else
                testCase.TestData.toTest{iTest, iArg} = testCase.TestData.(currentArg{1})(currentArg{2:end});
            end
        end
    end

    try
        % Run test.
        outputs = cell(1, nargout(testCase.TestData.functionToTest));
        [outputs{:}] = testCase.TestData.functionToTest(testCase.TestData.toTest{iTest, 1:end-2});

        if testCase.TestData.toTest{iTest, end-1}  % an error was expected but wasn't raised.
            result = 'failed';
            message = sprintf("error '%s' wasn't raised", testCase.TestData.toTest{iTest, end-1});
        else  % evaluate output
            try
                [result, message] = testCase.TestData.evaluateOutput(testCase.TestData.toTest{iTest, 1:end-2}, ...
                                                                     outputs, ...
                                                                     testCase.TestData.toTest{iTest, end});
            catch ME
                result = 'failed';
                message = sprintf('evaluation failed: %s - %s', ME.identifier, ME.message);
            end
        end

    catch ME  % the tested function raised an error
        if testCase.TestData.toTest{iTest, end-1}  % an error was expected
            if strcmp('error', testCase.TestData.toTest{iTest, end-1})
                result = 'passed';
                message = 'generic error detection';
            elseif strcmp(ME.identifier, testCase.TestData.toTest{iTest, end-1})
                result = 'passed';
                message = '';
            else
                result = 'warning';
                message = sprintf("raised '%s', but should have raised '%s'", ...
                                  ME.identifier, testCase.TestData.toTest{iTest, end-1});
            end
        else  % no error was expected
            result = 'failed';
            message = sprintf('unexpected error: %s: %s', ME.identifier, ME.message);
        end
    end

    % If warning or failed, let the user know.
    if ~strcmp(result, 'passed')
        fprintf('test %d %s: %s\n', iTest, result, message)
        pause(0.1);
    end
    
    % Debug pause.
    if isfield(testCase.TestData, 'debug') && testCase.TestData.debug > 0
        % Tips: The test failed, try to figure out why ->
        %       - look at the current test in testCase.TestData.toTest, test number: iTest
        %       - the output of the current test is in: outputs
        %       - run the function again:
        %           testCase.TestData.functionToTest(testCase.TestData.toTest{iTest, 1:end-2});
        %       - run the evaluation again:
        %           testCase.TestData.evaluateOutput(testCase.TestData.toTest{iTest, 1:end-2}, outputs, testCase.TestData.toTest{iTest, end});
        try
            if testCase.TestData.debug == 1 && strcmp(result, 'failed')
                error('EMC_runTest:debug', 'test failed')
            elseif testCase.TestData.debug == 2 && ~strcmp(result, 'passed')
                 error('EMC_runTest:debug', 'test failed')
            end
        catch
            continue
        end
    end

    % Save the results.
    results(iTest, :) = {result, message};

    % Clean fixture.
    if any(~cellfun(@isempty, cleanArg))
        for iArg = 1:length(cleanArg)
            if ~isempty(cleanArg{iArg})
                testCase.TestData.toTest{iTest, iArg} = cleanArg{iArg};
            end
        end
    end
end  % for every test

% Every tests were ran. Update the log with the results.
testCase.TestData.toTest = [testCase.TestData.toTest, results];
save(logFileName, 'testCase');

end  % EMC_runTest

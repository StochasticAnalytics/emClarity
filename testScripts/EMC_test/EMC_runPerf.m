function EMC_runPerf(testCase)
%
% Run perfomance tests (via timeit, gputimeit).
%
% Inputs:
% testCase.TestData should contains the following field:
%
%   functionToTest (handle):    Function to test.
%                               Signature: [ouput1, ...] = example(input1, ...)
%                               Functions without outputs and/or inputs are also accepted.
%
%   toTest (cell):            	Should correspond to functionToTest inputs (same order), plus 1 argument:
%                               - method (str): 'cpu' or 'gpu'.
%
%                               dimension:      cell(nTests, (nInputs + 1))
%                               format:         {input1A, input1B, input1C, ..., method;
%                                                input2A, input2B, input2C, ..., method;
%                                                ...
%                                                inputNA, inputNB, inputNC, ..., method}
%
%   (optional)
%   fixture_* (handle):         Same as EMC_runTest.
%
%   (optional)
%   testName (string):          Same as EMC_runTest.
%
% Notes:
%   -
%
% Examples:
%   One perfomance test with EMC_resize:
%   testCase.TestData.functionToTest = @EMC_resize;
%
%   % Create inputs with one fixtures.
%   method = 'gpu';  % test with gputimeit
%   testCase.TestData.fixtureImg = @(Size, Precision) ones(Size, Precision);
%   img = {'fixtureImg', method, 'single', [1280,1280]};  % the array will be created just before running the test.
%   limits = [10,-10,10,-10];
%   option = {};
%   
%   testCase.TestData.toTest = {img, limits, option, method};
%
%   % Run the test.
%   EMC_runPerf(testCase);  % logfile: ./logPerf/EMC_resize/EMC_resize.mat
%
% Other m-files required:
%
% Created: 18Jan2020
% Last revision: 18Jan2020
% See also EMC_runTest
%

%% Save inputs as log.
stack = dbstack;
if length(stack) == 1  % EMC_runPerf is called from the cmd line.
    error('This function should be called from a script, not from the command line.')
end
nameParentFunction = stack(2).name;
nameFunctionToTest = ['log_', func2str(testCase.TestData.functionToTest)];

folderName = [pwd, '/logPerf/', nameFunctionToTest];
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

%% Performance test
results = cell(nbTests, 1);  % execution times
for iTest = 1:nbTests  % for every test
    
    % Extract fixture.
    nbArgs = length(testCase.TestData.toTest(iTest, :)) - 1;  % minus method
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
    
    % Run performation test.
    if isstring(testCase.TestData.toTest{iTest, end}) || ischar(testCase.TestData.toTest{iTest, end})
        if strcmpi(testCase.TestData.toTest{iTest, end}, 'cpu')
            results{iTest, 1} = timeit(@()testCase.TestData.functionToTest(testCase.TestData.toTest{iTest, 1:end-1}));
        elseif strcmpi(testCase.TestData.toTest{iTest, end}, 'gpu')
            results{iTest, 1} = gputimeit(@()testCase.TestData.functionToTest(testCase.TestData.toTest{iTest, 1:end-1}));
        else
            results{iTest, 1} = "method should be 'cpu', 'cpu'";
        end
    else
        results{iTest, 1} = "method should be 'cpu', 'cpu'";
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
testCase.TestData.toTest = [testCase.TestData.toTest, results];
save(logFileName, 'testCase');

end  % EMC_runPerf

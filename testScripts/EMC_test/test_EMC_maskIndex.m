function tests = test_EMC_maskIndex
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_maskIndex;
testCase.TestData.evaluateOutput = @evaluateOutput;
end

function [result, message] = evaluateOutput(TYPE, SIZE, METHOD, OPTION, OUTPUTCELL, ~)
% check the size, method and precision
% if TYPE='fftshift', compare shift with fftshift
% if TYPE='ifftshift', compare shift with ifftshift
% if TYPE='nc2nc', compare with ifftn (same test as EMC_irfftn)

% if half SIZE, use the half2full option to recompute the full grid,
% then compare with fftshift or ifftshift.

result = 'failed';

if length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    INDEX = OUTPUTCELL{1};
end

% precision and method
[actualMethod, actualPrecision] = help_getClass(INDEX);

if ~strcmp(actualMethod, METHOD)
    message = sprintf('expected method=%s, got %s', METHOD, actualMethod);
    return
end

if help_isOptionDefined(OPTION, 'precision')
    expectPrecision = help_getOptionParam(OPTION, 'precision');
    if strcmpi(expectPrecision, 'uint')
        if prod(SIZE) <= 2^16
            expectPrecision = 'uint16';
        elseif prod(SIZE) <= 2^32
            expectPrecision = 'uint32';
        else
            expectPrecision = 'uint64';
        end
    elseif strcmpi(expectPrecision, 'int')
        if prod(SIZE) <= 2^16
            expectPrecision = 'int16';
        elseif prod(SIZE) <= 2^32
            expectPrecision = 'int32';
        else
            expectPrecision = 'int64';
        end
    end
else
    if prod(SIZE) <= 2^16
        expectPrecision = 'uint16';  % default
    elseif prod(SIZE) <= 2^32
        expectPrecision = 'uint32';  % default
    else
        expectPrecision = 'uint64';  % default
    end
end
if ~strcmp(actualPrecision, expectPrecision)
    message = sprintf('expected precision=%s, got %s', expectPrecision, actualPrecision);
    return
end

% size
if help_isOptionDefined(OPTION, 'half') && help_getOptionParam(OPTION, 'half')
    halfSize = [floor(SIZE(1)/2)+1, SIZE(2:end)];
    if ~isequal(halfSize, size(INDEX))
        message = sprintf('expected half size=%d, got %d', mat2str(halfSize), mat2str(size(INDEX)));
        return
    end
elseif size(INDEX) ~= SIZE
 	message = sprintf('expected size=%d, got %d', mat2str(SIZE), mat2str(size(INDEX)));
   	return
end


if strcmpi(TYPE, 'fftshift')
    if help_isOptionDefined(OPTION, 'half') && help_getOptionParam(OPTION, 'half')
        % shift back and test if you have the same img.
        testImg = rand(halfSize);
        centeredImg = testImg(INDEX);
        non_centeredImg = centeredImg(EMC_maskIndex('ifftshift', SIZE, 'cpu', {'half', true}));
        if any(abs(testImg - non_centeredImg) > 1e-10, 'all')
            message = 'fftshift is not reverse by ifftshift';
            return
        end
    else  % full grid; compare with fftshift
        testImg = rand(SIZE);
        if any(abs(testImg(INDEX) - fftshift(testImg)) > 1e-10, 'all')
            message = 'shift is not equal to fftshift';
            return
        end
    end
elseif strcmpi(TYPE, 'ifftshift')
    if help_isOptionDefined(OPTION, 'half') && help_getOptionParam(OPTION, 'half')
        % shift back and test if you have the same img.
        testImg = rand(halfSize);
        non_centeredImg = testImg(INDEX);
        centeredImg = non_centeredImg(EMC_maskIndex('fftshift', SIZE, 'cpu', {'half', true}));
        if any(abs(testImg - centeredImg) > 1e-10, 'all')
            message = 'fftshift is not reverse by ifftshift';
            return
        end
    else  % full grid; compare with ifftshift
        testImg = rand(SIZE);
        if any(abs(testImg(INDEX) - ifftshift(testImg)) > 1e-10, 'all')
            message = 'shift is not equal to ifftshift';
            return
        end
    end
end

% nc2nc is tested in test_EMC_rfftn_and_irfftn

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default(testCase)
types = {'fftshift'; 'ifftshift'; 'nc2nc'; "fftshift"; "ifftshift"; "nc2nc"};
sizes = [help_getRandomSizes(5, [500, 2000], '2d'); help_getRandomSizes(5, [50, 200], '3d')];
methods = {'cpu'; 'gpu'};
options = help_getBatchOption({'half', {true; false}; ...
                               'precision', {'int'; 'uint'; 'single'; 'double'}});
testCase.TestData.toTest = help_getBatch(types, sizes, methods, options, {false}, {false});
EMC_runTest(testCase);

end


function test_assumptions(testCase)

% type
types = {'kayak'; ''; nan; inf; [1,2,3]; {''};};
testCase.TestData.toTest = help_getBatch(types, {[10,10]}, {'cpu'}, {{}}, {'EMC:TYPE'}, {false});

% size
sizes = {'kayak'; [1,10]; [10,1]; [10,10,10,10]; []; 2; {}; [0,10]; [10.5;10]; ...
         [-10,10]; [nan, 10]; [inf,10]; [10;10]};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
	help_getBatch({'fftshift'}, sizes, {'cpu'}, {{}}, {'EMC:SIZE'}, {false})];

% method
methods = {'kayak'; 1; {}; [1,2]; []; nan; inf};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'fftshift'}, {[10,10]}, methods, {{}}, {'EMC:METHOD'}, {false})];

% option
options = help_getBatchOption({'precision', {'kayak'; 12; []; nan; inf}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'fftshift'}, {[10,10]}, {'cpu'}, options, {'EMC:precision'}, {false})];

options = help_getBatchOption({'half', {'kayak'; 12; []; nan; inf; 1; 0}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
    help_getBatch({'fftshift'}, {[10,10]}, {'cpu'}, options, {'EMC:half'}, {false})];

EMC_runTest(testCase);
end

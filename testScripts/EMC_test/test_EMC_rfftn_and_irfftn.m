function tests = test_EMC_rfftn_and_irfftn
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 0;
testCase.TestData.functionToTest = @EMC_irfftn_and_rfftn;
testCase.TestData.evaluateOutput = @evaluateOutput;
testCase.TestData.fixtureRand = @help_getInputRand;
end

function [dft, img] = EMC_irfftn_and_rfftn(IMAGE)
dft = EMC_rfftn(IMAGE);
img = EMC_irfftn(dft, size(IMAGE));
end


function [result, message] = evaluateOutput(IMAGE, OUTPUTCELL, ~)
% Check that precision and method match the IMAGE.
% Check dft sX is half.
% Check EMC_irfftn(EMC_rfftn(img), size(img)) is equal to ifftn(fftn(img)).
%

result = 'failed';

if length(OUTPUTCELL) ~= 2
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    [dft, outImg] = OUTPUTCELL{:};
end

% precision and method
[actualMethod1, actualPrecision1] = help_getClass(dft);
[actualMethod2, actualPrecision2] = help_getClass(outImg);
[expectMethod, expectPrecision] = help_getClass(IMAGE);

if ~strcmp(actualMethod1, expectMethod)
    message = sprintf('expected method=%s, got %s', expectMethod, actualMethod1);
    return
elseif ~strcmp(actualPrecision1, expectPrecision)
    message = sprintf('expected precision=%s, got %s', expectPrecision, actualPrecision1);
    return
elseif ~strcmp(actualMethod2, expectMethod)
    message = sprintf('expected method=%s, got %s', expectMethod, actualMethod2);
    return
elseif ~strcmp(actualPrecision2, expectPrecision)
    message = sprintf('expected precision=%s, got %s', expectPrecision, actualPrecision2);
    return
end

% half size
dftSize = size(dft);
inSize = size(IMAGE);
halfSize = [floor(inSize(1)/2)+1, inSize(2:end)];
if ~isequal(halfSize, dftSize)
    message = sprintf('expected size=%s, got %s', ...
                      mat2str(halfSize), mat2str(dftSize));
    return
end

% size
if ~isequal(size(outImg), inSize)
    message = sprintf('expected size=%s, got %s', mat2str(size(outImg)), mat2str(inSize));
    return
end

% output image should be equal to input image.
if any(abs(outImg - IMAGE) > 1e-5, 'all')
 	message = 'output image is not equal to input image';
	return
end

result = 'passed';
message = '';

end


function test_default(testCase)
sizes = [help_getRandomSizes(50, [10,300], '2d'); help_getRandomSizes(50, [10,200], '3d')];
img = help_getBatch({'fixtureRand'}, {'cpu';'gpu'}, {'double';'single'}, sizes);
testCase.TestData.toTest = help_getBatch(img, {false}, {false});

EMC_runTest(testCase);

end

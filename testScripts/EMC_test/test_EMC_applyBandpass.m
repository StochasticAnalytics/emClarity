function tests = test_EMC_applyBandpass
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @EMC_applyBandpass;
testCase.TestData.evaluateOutput = @evaluateOutput;
testCase.TestData.fixtureRandn = @getRandomImgMean;
end


function img = getRandomImgMean(method, precision, size)
img = help_getInputRand(method, precision, size) .* (rand(1) * 10);
img = img ./ max(img(:));
end


function [img1, img2] = compareFullAndHalf(img, ~, option)

[METHOD, PRECISION] = help_getClass(img);
bp_half = EMC_getBandpass(size(img), 1, 15, 2.5, METHOD, {'precision', PRECISION; 'half', true});
bp_full = EMC_getBandpass(size(img), 1, 15, 2.5, METHOD, {'precision', PRECISION; 'half', false});

img1 = EMC_applyBandpass(img, bp_half, option);
img2 = EMC_applyBandpass(img, bp_full, option);

% img1 and img2 should be identical.
end


function [result, message] = evaluateOutput(IMAGE, ~, OPTION, OUTPUTCELL, ~)
% IMAGE = EMC_applyBandpass(IMAGE, BANDPASS, OPTION)
% Check:
%   METHOD:     output has the same method has IMAGE
%   PRECISION:  output has the same precision has IMAGE
%   SIZE:       output has the same size has IMAGE
%   Output should be real if 'iff'=true and complex if false.
%   Mean and std should be close to 0 and 1 respectively, if 'standardize' is true.
%
%   If EXTRA, load and check for equality with fixture.
%

result = 'failed';

if length(OUTPUTCELL) == 2
    [filteredImg, filteredImg2] = OUTPUTCELL{:};
    ishalf = true;
elseif length(OUTPUTCELL) ~= 1
    message = sprintf('expected output number = 1, got %d', length(OUTPUTCELL));
    return
else
    filteredImg = OUTPUTCELL{1};
    ishalf = false;
end

% precision and method
[actualMethod, actualPrecision] = help_getClass(filteredImg);
[expectMethod, expectPrecision] = help_getClass(IMAGE);

if ~strcmp(actualMethod, expectMethod)
    message = sprintf('expected method=%s, got %s', expectMethod, actualMethod);
    return
elseif ~strcmp(actualPrecision, expectPrecision)
    message = sprintf('expected precision=%s, got %s', expectPrecision, actualPrecision);
    return
end

% ifft
if help_isOptionDefined(OPTION, 'ifft')
    if help_getOptionParam(OPTION, 'ifft')
        isifft = true;
        if ~isreal(filteredImg)
            message = 'output should be real with because ifft=true';
            return
        % size
        elseif size(filteredImg) ~= size(IMAGE)
            message = sprintf('expected size=%s, got %s', mat2str(size(IMAGE)), mat2str(size(filteredImg)));
            return
        end
    else
        isifft = false;
        if isreal(filteredImg)
            message = 'output should be complex with because ifft=true';
            return
        % size
        elseif ishalf
            imgSize = size(IMAGE);
            if ~isequal([floor(imgSize(1)/2)+1, imgSize(2:end)], size(filteredImg))  % half size
                message = sprintf('expected size=%s, got %s', ...
                    mat2str([floor(imgSize(1)/2)+1, imgSize(2:end)]), mat2str(size(filteredImg)));
                return
            end
        elseif size(filteredImg) ~= size(IMAGE)
            message = sprintf('expected size=%s, got %s', mat2str(size(IMAGE)), mat2str(size(filteredImg)));
            return
        end 
    end
else
    isifft = true;
    if ~isreal(filteredImg)
        message = 'output should be real with because ifft=true';
       	return
   % size
    elseif size(filteredImg) ~= size(IMAGE)
        message = sprintf('expected size=%s, got %s', mat2str(size(IMAGE)), mat2str(size(filteredImg)));
        return
    end
end

% uniform
if isreal(filteredImg)
    if help_isOptionDefined(OPTION, 'standardize')
        if help_getOptionParam(OPTION, 'standardize')
            if abs(mean(filteredImg(:))) > 0.1
                message = sprintf('output mean should be 0, got %f, original %f', mean(filteredImg(:)), mean(IMAGE(:)));
                return
            elseif abs(std(filteredImg(:)) - 1) > 0.1
                message = sprintf('output std should be 1, got %f', std(filteredImg(:)));
                return
            end
        end
    else
        if abs(mean(filteredImg(:))) > 0.1
            message = sprintf('output mean should be 0, got %f, original %f', mean(filteredImg(:)), mean(IMAGE(:)));
            return
        elseif abs(std(filteredImg(:)) - 1) > 0.1
           message = sprintf('output std should be 1, got %f', std(filteredImg(:)));
           return
        end
    end
end

% compare with fixture
if ishalf && isifft
    if any(abs(filteredImg - filteredImg2) > 1e-5, 'all')
       	message = 'half not equivalent to full';
       	return
    end
end

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default_2d(testCase)
% Test for 2d, cpu/gpu and every option.
%

sizes = help_getRandomSizes(2, [500, 1000], '2d');
precisions = {'single'; 'double'};
methods = {'cpu'; 'gpu'};
img = help_getBatch({'fixtureRandn'}, methods, precisions, sizes);

option = help_getBatchOption({'ifft', {true;false}; ...
                              'standardize', {true;false}});

testCase.TestData.toTest = help_getBatch(img, option, {false}, {false});

% bandpass should have same precision and method than image
toTest = cell(length(testCase.TestData.toTest), 5);
toTest(:,1) = testCase.TestData.toTest(:,1);
toTest(:,2) = testCase.TestData.toTest(:,1);
toTest(:,3:end) = testCase.TestData.toTest(:,2:end);

testCase.TestData.toTest = toTest;

EMC_runTest(testCase);

end


function test_default_3d(testCase)
% Test for 3d, cpu/gpu and every option.
%

sizes = help_getRandomSizes(2, [150,250], '3d');
precisions = {'single'; 'double'};
methods = {'cpu'; 'gpu'};
img = help_getBatch({'fixtureRandn'}, methods, precisions, sizes);

option = help_getBatchOption({'ifft', {true;false}; ...
                              'standardize', {true;false}});

testCase.TestData.toTest = help_getBatch(img, option, {false}, {false});

% bandpass should have same precision and method than image
toTest = cell(length(testCase.TestData.toTest), 5);
toTest(:,1) = testCase.TestData.toTest(:,1);
toTest(:,2) = testCase.TestData.toTest(:,1);
toTest(:,3:end) = testCase.TestData.toTest(:,2:end);

testCase.TestData.toTest = toTest;

EMC_runTest(testCase);

end


function test_assumptions(testCase)
% make sure error are raised

% image
imgIm = ones(10,10) + 1i;
img = {'kayak'; ones(5,5,5,5); []; 2; {}; imgIm};
testCase.TestData.toTest = help_getBatch(img, {ones(5,5)}, {{}}, {'error'}, {false});

% bandpass
bandpass = {'kayak'; ones(1,10); ones(10,1); ones(5,5,5,5); []; 2; {}};
testCase.TestData.toTest = [testCase.TestData.toTest; ...
                            help_getBatch({ones(10,10)}, bandpass, {{}}, {'error'}, {false})];

testCase.TestData.toTest = [testCase.TestData.toTest; ...
    {ones(10,10,'single'), ones(10,10,'double'),             {}, 'EMC:IMAGE', false; ...
     ones(10,10,'double'), ones(10,10,'single'),             {}, 'EMC:IMAGE', false; ...
     ones(10,10,'single', 'gpuArray'), ones(10,10,'double'), {}, 'EMC:IMAGE', false; ...
     ones(10,10,'double', 'gpuArray'), ones(10,10,'double'), {}, 'EMC:IMAGE', false; ...
     ones(10,10,'single'), ones(10,10,'double', 'gpuArray'), {}, 'EMC:IMAGE', false; ...
     ones(10,10,'double'), ones(10,10,'double', 'gpuArray'), {}, 'EMC:IMAGE', false; ...
     ones(10,10,'single', 'gpuArray'), ones(10,10,'double', 'gpuArray'), {}, 'EMC:IMAGE', false; ...
     ones(10,10,'double', 'gpuArray'), ones(10,10,'single', 'gpuArray'), {}, 'EMC:IMAGE', false; ...
     ones(10,10), ones(5,10), {}, 'EMC:IMAGE', false}];

% options
options = help_getBatchOption({'iff', {1;0;-1;3;'kayak';nan;inf;{};[];''}; ...
                               'standardize', {1;0;-1;3;'kayak';inf;nan;{};[];''}});
options = options(~cellfun(@isempty, options), :);
testCase.TestData.toTest = [testCase.TestData.toTest; ...
        help_getBatch({ones(10,10)}, ones(10,10), options, {'error'}, {false})];
    
% check that bandpass other than full and half raise an error.
testCase.TestData.toTest = [testCase.TestData.toTest; ...
       {ones(100,100), ones(80,100), {},    'EMC:IMAGE', false; ...
       ones(50,50,50), ones(25,50,50), {}, 'EMC:IMAGE', false}];

EMC_runTest(testCase);

end


function test_half(testCase)
% check that using half bandpass gives the same filtered image at the end than full bandpass.
testCase.TestData.functionToTest = @compareFullAndHalf;

sizes = help_getRandomSizes(2, [150,250], '2d');
precisions = {'single'; 'double'};
methods = {'cpu'; 'gpu'};
img = help_getBatch({'fixtureRandn'}, methods, precisions, sizes);
option = help_getBatchOption({'ifft', {true;false}; 'standardize', {true;false}});
testCase.TestData.toTest =  help_getBatch(img, {nan}, option, {false}, {false});

sizes = help_getRandomSizes(2, [150,250], '3d');
img = help_getBatch({'fixtureRandn'}, methods, precisions, sizes);
testCase.TestData.toTest = [testCase.TestData.toTest; help_getBatch(img, {nan}, option, {false}, {false})];

EMC_runTest(testCase);

end

function tests = test_EMC_maskReference
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.debug = 2;
testCase.TestData.functionToTest = @p_EMC_maskReference;
testCase.TestData.evaluateOutput = @evaluateOutput;
testCase.TestData.fixture_Rand = @help_getInputRand;

global EMC_gp_verbose
EMC_gp_verbose = false;
end


function [o1, o2, o3, o4] = p_EMC_maskReference(IMAGE, PIXEL, OPTION)

if help_isOptionDefined(OPTION, 'fsc') && help_getOptionParam(OPTION, 'fsc')
    [o1, o2, o3, o4] = EMC_maskReference(IMAGE, PIXEL, OPTION);
else
    [o1, o2] = EMC_maskReference(IMAGE, PIXEL, OPTION);
    o3 = 'empty';
    o4 = 'empty';
end

end


function [result, message] = evaluateOutput(IMAGE, ~, OPTION, OUTPUTCELL, ~)
% [MASK, COM]                      = EMC_maskReference(IMAGE, PIXEL, OPTION) if 'fsc'=false.
% [MASK, MASK_CORE, FRACTION, COM] = EMC_maskReference(IMAGE, PIXEL, OPTION) if 'fsc'=true.
%
% Check:
%   METHOD:     MASK, MASK_CORE, COM are following the METHOD instruction.
%   PRECISION:  MASK, MASK_CORE, COM have the desired precision
%   SIZE:       MASK, MASK_CORE, COM, FRACTION have the correct size.
%
%   If EXTRA, load and check for equality with fixture.
%

result = 'failed';

if help_isOptionDefined(OPTION, 'fsc') && help_getOptionParam(OPTION, 'fsc')
    isFsc = true;
else
    isFsc = false;
end

% method and precision
[expectedMethod, expectedPrecision] = help_getClass(IMAGE);

% MASK - precision, method and size
[method, precision] = help_getClass(OUTPUTCELL{1});
if ~strcmpi(method, expectedMethod)
    message = sprintf('MASK, expected method=%s, got %s', expectedMethod, method);
    return
elseif ~strcmpi(precision, expectedPrecision)
    message = sprintf('MASK, expected precision=%s, got %s', expectedPrecision, precision);
    return
elseif ~isequal(size(IMAGE), size(OUTPUTCELL{1}))
    message = sprintf('MASK, expected size=%s, got %s',mat2str(size(IMAGE)), mat2str(size(OUTPUTCELL{1})));
    return
end

if isFsc
    % MASK_CORE - precision, method and size
    [method, precision] = help_getClass(OUTPUTCELL{2});
    if ~strcmp(method, expectedMethod)
        message = sprintf('MASK_CORE, expected method=%s, got %s', expectedMethod, method);
        return
    elseif ~strcmp(precision, expectedPrecision)
        message = sprintf('MASK_CORE, expected precision=%s, got %s', expectedPrecision, precision);
        return
    elseif ~isequal(size(IMAGE), size(OUTPUTCELL{2}))
        message = sprintf('MASK_CORE, expected size=%s, got %s', mat2str(size(IMAGE)), mat2str(size(OUTPUTCELL{2})));
        return
    end
    
    % FRACTION
    if ~isnumeric(OUTPUTCELL{3}) || ~isscalar(OUTPUTCELL{3}) || isinf(OUTPUTCELL{3}) || isnan(OUTPUTCELL{3})
        message = sprintf('FRACTION should be a numeric scalar, not inf and not nan');
        return
    end
    
    % COM - precision, method and size
    if help_isOptionDefined(OPTION, 'com') && help_getOptionParam(OPTION, 'com')
        [method, precision] = help_getClass(OUTPUTCELL{4});
        if ~strcmp(method, expectedMethod)
            message = sprintf('COM, expected method=%s, got %s', expectedMethod, method);
            return
        elseif ~strcmp(precision, expectedPrecision)
            message = sprintf('COM, expected precision=%s, got %s', expectedPrecision, precision);
            return
        elseif ~isequal(size(size(IMAGE)), size(OUTPUTCELL{4}))
            message = sprintf('COM, expected size=%s, got %s', mat2str(size(size(IMAGE))), mat2str(size(OUTPUTCELL{4})));
            return
        end
    else
        if ~isscalar(OUTPUTCELL{4}) || ~isnan(OUTPUTCELL{4})
            message = sprintf('COM should be nan');
            return
        end
    end
else
    % COM - precision, method and size
    if help_isOptionDefined(OPTION, 'com') && help_getOptionParam(OPTION, 'com')
        [method, precision] = help_getClass(OUTPUTCELL{2});
        if ~strcmp(method, expectedMethod)
            message = sprintf('COM, expected method=%s, got %s', expectedMethod, method);
            return
        elseif ~strcmp(precision, expectedPrecision)
            message = sprintf('COM, expected precision=%s, got %s', expectedPrecision, precision);
            return
        elseif ~isequal(size(size(IMAGE)), size(OUTPUTCELL{2}))
            message = sprintf('COM, expected size=%s, got %s', mat2str(size(size(IMAGE))), mat2str(size(OUTPUTCELL{2})));
            return
        end
    else
        if ~isscalar(OUTPUTCELL{2}) || ~isnan(OUTPUTCELL{2})
            message = sprintf('COM should be nan');
            return
        end
    end
end

% You have passed the test, congratulation.
result = 'passed';
message = '';

end


function test_default(testCase)
% Test for 2d/3d, cpu/gpu and every option.
%

sizes = [help_getRandomSizes(1, [50, 100], '2d'); help_getRandomSizes(1, [50, 100], '3d')];
precisions = {'single'; 'double'};
methods = {'cpu'; 'gpu'};
imgs = help_getBatch({'fixture_Rand'}, methods, precisions, sizes);

pixel = {1};

option = help_getBatchOption({'precision', {'single' ; 'double'}; ...
                              'method', {'cpu'; 'gpu'}; ...
                              'com', {true; false}; ...
                              'fsc', {true; false}; ...
                              });
option = [option; ...
          help_getBatchOption({'origin', {0; 1; 2}; ...
                               'hydration_scaling', {2.5}; ...
                               'lowpass', {20}; ...
                               'threshold', {3}; ...
                               'com', {true}; ...
                               'fsc', {true}; ...
                               })];

testCase.TestData.toTest = help_getBatch(imgs, pixel, option, {false}, {false});
EMC_runTest(testCase);

end

function BATCH = help_getBatch(varargin)
%
%
%
%

numberOfRuns = 1;
for iArg = 1:nargin
    [row, column] = size(varargin{iArg});
    numberOfRuns = numberOfRuns * row;
    
    % if a cell has more than one column, concatenate these columns into one 1xn cell.
    if column > 1
        tmp = cell(row, 1);
        for i=1:row; tmp(i, 1) = {varargin{iArg}(i, :)}; end
        varargin{iArg} = tmp;
    end
end

BATCH = cell(numberOfRuns, nargin);

for iArg = 1:nargin
    
    nargs = length(varargin{iArg});
    if iArg > 1
        repeat = repeat / nargs;
    else
        repeat = numberOfRuns / nargs;
    end
    loop = numberOfRuns / (repeat * nargs);
    
    count = 1;
    for iLoop = 1:loop    
        for iItem = 1:nargs
            for iRepeat = 1:repeat
                BATCH(count, iArg) = varargin{iArg}(iItem);
                count = count + 1;
            end
        end
    end
end

end
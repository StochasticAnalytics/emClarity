function BATCH = help_getBatchOption(OPTIONS)
%
% Compute every combinaison of optional parameters.
%

nbArgs = length(OPTIONS(:, 1));
inGetBatch = cell(1, nbArgs);
for iOp = 1:length(OPTIONS(:, 1))
    inGetBatch{1, iOp} = [{{nan}}; OPTIONS{iOp, 2}];
end

out = help_getBatch(inGetBatch{:});

% concatenate row after removing {nan}
BATCH = cell(length(out), 1);
for iOp = 1:length(out)
    row = out(iOp, :);
    
    % get indexes of {nan}
    valideIdx = true(1, numel(row));
    for i = 1:numel(row)
        if iscell(row{i}) && ~isempty(row{i}) && length(row{i}) == 1 && any(isnan(row{i}{1}), 'all')
            valideIdx(i) = 0;
        end
    end
    
    % create final option cell
    if sum(valideIdx) > 0
        intermediateOut = cell(sum(valideIdx), 2);
        count = 1;
        for iParam = 1:numel(valideIdx)
            if valideIdx(iParam) == 1
                intermediateOut{count, 1} = OPTIONS{iParam, 1};
                intermediateOut{count, 2} = row{1, iParam};
                count = count + 1;
            end
        end
        BATCH{iOp, 1} = intermediateOut;
    else
        BATCH{iOp, 1} = {};
    end
end

end
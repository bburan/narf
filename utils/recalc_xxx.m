function recalc_xxx(start_depth)
% Recalculates the XXX data structure from start_depth onwards
global XXX STACK;

% If data already exists, invalidate it just in case
if start_depth < length(XXX)
    XXX = XXX(1:start_depth);  % Invalidate later data so it cannot be used
end

% If there is not enough data, begin the calculation from the top
if start_depth > length(XXX)
    start_depth = length(XXX);
end

% Now, do the recalculation of the data
for ii = start_depth:length(STACK);
    if ~STACK{ii}.isready_pred(STACK(1:ii), XXX(1:ii));
        error('Stack was not fully ready at depth %d', ii);
    end
    XXX{ii+1} = STACK{ii}.fn(STACK(1:ii), XXX(1:ii));
end

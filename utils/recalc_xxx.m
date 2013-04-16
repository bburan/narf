function recalc_xxx(start_depth, end_depth)
% recalc_xxx(start_depth, end_depth)
%
% Recalculates the XXX data structure from = start_depth to < end_depth,
% invalidating and overwriting any previous data in that area.
% A very important function, recalc_xxx is responsible for actually 
% computing the prediction made by a model.
%
% ARGUMENTS:
%    start_depth    The point of the stack to begin recomputation from.
%                   Must be in range 1 <= start_depth <= length(STACK)
%    end_depth      Defaults to the end of the stack.
% 
% RETURNS: Nothing
% 

global XXX STACK;

if start_depth < 1
    error('Start depth less than 1 not allowed.');
end

if ~exist('end_depth', 'var')
    end_depth = length(STACK);
end

% If data already exists, invalidate it just in case
if start_depth < length(XXX)
    XXX = XXX(1:start_depth);  % Invalidate later data so it cannot be used
end

% If there is not enough data, begin the calculation from the top
if start_depth > length(XXX)
    start_depth = length(XXX);
end

% Now, do the recalculation of the data
for ii = start_depth:(end_depth-1),
    m = STACK{ii};
    
    if iscell(m)        
        % It is a split parameter module. Split the data for each parameter set
        splitxxx = m{1}.splitter(XXX(1:ii));        
        
        ret = cell(1, length(m));
        
        % Check sizes for consistency
        if length(splitxxx) ~= length(m)
            error('Splitter group count does not match STACK parameter set count.');
        end
        
        % Iterate through each parameter set
        for jj = 1:length(m)
            tmpstack = STACK(1:ii-1);
            tmpstack{end+1} = m{jj};         

            if ~m{jj}.isready_pred(tmpstack, splitXXX(jj));
                error('Stack was not fully ready at depth %d idx %d', ii, jj);
            end
            ret{end+1} = m{jj}.fn(m, splitxxx{jj}{ii}, tmpstack, splitxxx(jj));
        end
        
        % Unify the returned values 
        XXX{ii+1} = m{1}.unifier(ret);
    else
        % It is not a split parameter module
        if ~STACK{ii}.isready_pred(STACK(1:ii), XXX(1:ii));
            error('Stack was not fully ready at depth %d', ii);
        end
        XXX{ii+1} = m.fn(m, XXX{ii}, STACK(1:ii), XXX(1:ii));
    end
end

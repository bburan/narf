function calc_xxx(start_depth, end_depth)
% calc_xxx(start_depth, end_depth)
%
% Computes the XXX data structure from = start_depth to < end_depth,
% invalidating and overwriting any previous data in that area.
% A very important function, calc_xxx is responsible for actually 
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

% Compute the entire stack, using calc_paramsets to manage actual
% computation of the modules' parameter sets
for ii = start_depth:(end_depth-1),
    [mdls, ~, xouts] = calc_paramsets(STACK(1:ii), XXX(1:ii));
    
    if length(xouts) == 1 && ~isfield(mdls{1}, 'unifier')
        XXX{ii+1} = xouts{1};
    elseif isfield(mdls{1}, 'unifier')
        XXX{ii+1} = mdls{1}.unifier(xouts);
    else
        error('A unifier should not be used for single parameter set.');
    end    
end

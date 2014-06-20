function calc_xxx(start_depth, end_depth)
% calc_xxx(start_depth)
%
% Computes the XXX data structure, overwriting any previous data as needed.
%
% Now a wrapper around update_xxx(). 
%
% ARGUMENTS:
%    start_depth    The depth in the stack to begin recomputation from.
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
%if start_depth < length(XXX)
%    XXX = XXX(1:start_depth);  % Invalidate later data so it cannot be used
%end

% If there is not enough data, begin the calculation from the top
if start_depth > length(XXX)
    start_depth = length(XXX);
end

% If the end depth is beyond the length of the stack, stop
if end_depth > length(STACK)
    end_depth = length(STACK);
end

flat_start_depth = 1;
flat_end_depth = 1;
for ii = 1:length(STACK)
    for jj = 1:length(STACK{ii})
        if ii < start_depth
            flat_start_depth = flat_start_depth + 1;
        end
        if ii <= end_depth
            flat_end_depth = flat_end_depth + 1;
        end
    end
end

update_xxx(flat_start_depth, flat_end_depth - 1);
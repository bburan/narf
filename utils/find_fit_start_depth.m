function [start_depth, end_depth]= find_fit_start_depth(stack)
% [start_depth, end_depth] = find_fit_start_depth(stack)
%
% Returns the indexes of STACK containing first & last fittable fields. 
% This is important to reduce the amount of work that fitting algorithms do, since
% they only need to recalculate the stack from the point where fields are
% actually fit. Very rarely, you may also need to know where the fittable
% fields end (such as during jackknifing), so results can be merged across
% jackknifes and passed through a performance metric.

start_depth = NaN;
end_depth = NaN;

for d = 1:length(stack)
    m = stack{d}{1};
    
    if isfield(m, 'fit_fields') && ~isempty(m.fit_fields)
        if isnan(start_depth)
            start_depth = d;
        end
        end_depth = d;
    end
end

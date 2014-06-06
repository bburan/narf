function [start_depth, end_depth]= find_fit_start_depth(~)
% [start_depth, end_depth] = find_fit_start_depth()
%
% Returns the module indexes (NOT the same as the STACK index) that
% contain the first & last fittable fields. 
% This is important to reduce the work that fitting algorithms do, since
% they only need to recalculate the stack from the point where fields are
% actually fit. Very rarely, you may also need to know where the fittable
% fields end (such as during jackknifing), so results can be merged across
% jackknifes and passed through a performance metric.

global STACK;

start_depth = NaN;
end_depth = NaN;

depth = 1;
for ii = 1:length(STACK)
    mm = STACK{ii};
    nsplits = length(mm);
    for kk = 1:nsplits
        m = mm{kk};                
        if isnan(start_depth) && isfield(m, 'fit_fields') && ~isempty(m.fit_fields)
            start_depth = depth;
        end
        end_depth = depth;
    end
    depth = depth + 1;
end
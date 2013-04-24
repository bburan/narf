function d = find_fit_start_depth(stack)
% d = find_fit_start_depth(stack)
%
% Returns the index of STACK containing the first fittable field. This is
% important to reduce the amount of work that fitting algorithms do, since
% they only need to recalculate the stack from the point where fields are
% actually fit. 
% 
for d = 1:length(stack)
    m = stack{d}{1};
    
    if isfield(m, 'fit_fields') && ~isempty(m.fit_fields)
        return;
    end
end

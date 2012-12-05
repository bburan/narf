function d = find_fit_start_depth(stack)
% Find the depth at which to start recalculating the stack
for d = 1:length(stack)
    if isfield(stack{d}, 'fit_fields') && ~isempty(stack{d}.fit_fields)
        return;
    end
end

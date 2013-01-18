function [mod, mod_idx] = find_module(stack, mod_name)
% Return the module with name 'mod_name'. 

for idx = 1:length(stack)
    if isequal(stack{idx}.name, mod_name)
        mod = stack{idx};
        mod_idx = idx;
        return
    end
end
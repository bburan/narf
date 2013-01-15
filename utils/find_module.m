function mod = find_module(stack, mod_name)
% Return the module with name 'mod_name'. 

for idx = 1:length(stack)
    if isequal(stack{idx}.name, mod_name)
        mod = stack{idx};
        return
    end
end
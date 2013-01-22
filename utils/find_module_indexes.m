function indexes = find_module_indexes(stack, mod_name)
% Return the indexes of all modules named mod_name 

indexes = [];
for idx = 1:length(stack)
    if isequal(stack{idx}.name, mod_name)
        indexes(end+1) = idx;
    end
end
function [mods, mod_idxs] = find_modules(stack, mod_name, first_only)
% [mods, mod_idxs] = find_modules(stack, mod_name, first_only=false)
%
% Returns a cell array of modules named mod_name, as well as their
% positional indexes in STACK. If first_only is provided, instead of a cell
% array being returned, the first matching module and index are returned.
%
% ARGUMENTS:
%    stack         The stack of modules
%    mod_name      A the module's string name
%
% RETURNS:
%    mods          A cell array of modules whose names matched mod_name
%    mod_idxs      A cell array of module indexes whose names matched
%
mods = {};
mod_idxs = {};

for idx = 1:length(stack)
    if isequal(stack{idx}.name, mod_name)
        if first_only
            mods = stack{idx};
            mod_idxs = idx;
            return;
        end
        mods{end+1} = stack{idx};
        mod_idxs{end+1} = idx;
    end
end

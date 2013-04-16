function append_module(module)
% APPEND_MODULE(module)
%
% Append module to the current global cell array STACK, then call the
% module's auto_init() method if it exists so that it automatically
% initializes itself. If the module has a 'splitter' and a 'unifier'
% function, then it actually appends a cell array.
%
% PLEASE use this instead of manually appending to STACK. Reasons why:
%   1. You don't need to care about module index numbers or parameter sets.
%   2. A module can auto-determine the size of an array (for example,
%   number of FIR coefficient channels) to match prior modules.

global STACK XXX;

l = length(STACK);

if all(isfield(module, {'splitter', 'unifier'}))    
    splitxxx = module.splitter(XXX(1:l+1));
    ret = {};
    
    for ii = 1:length(splitxxx)
        if isfield(module, 'auto_init')
            ret{end+1} = module.auto_init(STACK, splitxxx{ii});
        else
            ret{end+1} = module;
        end
    end
    
    STACK{l+1} = ret;
else
    if isfield(module, 'auto_init')
        STACK{l+1} = {module.auto_init(STACK, XXX(1:l))};
    else
        STACK{l+1} = {module};
    end
end

recalc_xxx(l+1); % was: XXX{l+2} = STACK{l+1}.fn(STACK(1:l+1), XXX(1:l+1));

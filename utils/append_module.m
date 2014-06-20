function append_module(module)
% APPEND_MODULE(module)
%
% Append module to the current global cell array STACK, then call the
% module's auto_init() method if it exists so that it automatically
% initializes itself. 
%
% PLEASE use this instead of manually appending to STACK. Reasons why:
%   1. You don't need to care about module index numbers or parameter sets.
%   2. A module can auto-determine the size of an array (for example,
%   number of FIR coefficient channels) to match prior modules.

global STACK XXX;

l = length(STACK);

if isfield(module, 'auto_init')
    STACK{l+1} = {module.auto_init(STACK, XXX(1:l+1))};
else
    STACK{l+1} = {module};
end

update_xxx(length(get_flatstack()));

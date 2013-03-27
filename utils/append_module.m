function append_module(module)
% APPEND_MODULE(module)
%
% Append module to the current global cell array STACK, then call the
% module's auto_init() method if it exists so that it automatically
% initializes itself. 
%
% Intended to be used when:
%   1. You don't care about module index numbers
%   2. A module needs to determine the size of an array (for example,
%   number of FIR coefficient channels) to match prior things.

global STACK XXX;

l = length(STACK);

if isfield(module, 'auto_init')
    STACK{l+1} = module.auto_init(STACK, XXX(1:l+1));
else
    STACK{l+1} = module;
end

XXX{l+2} = STACK{l+1}.fn(STACK(1:l+1), XXX(1:l+1));

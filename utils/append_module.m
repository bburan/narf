function append_module(module)
% Append module to the current STACK, and ask the module to use XXX data
% to initialize itself automatically. Good for script-writing when
% you don't want to keep track of the index number of a particular module,
% but want to build a model automatically.

global STACK XXX;

l = length(STACK);

if isfield(module, 'auto_init')
    STACK{l+1} = module.auto_init(STACK, XXX(1:l+1));
else
    STACK{l+1} = module;
end

XXX{l+2} = STACK{l+1}.fn(STACK(1:l+1), XXX(1:l+1));

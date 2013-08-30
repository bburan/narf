function scaat()
global STACK;
% Only append MSE module if needed
mods = find_modules(STACK, 'mean_squared_error', true);
if isempty(mods)
    mse();    
end
fit_scaat();
function nmse()

global MODULES META STACK;

mods = find_modules(STACK, 'mean_squared_error', true);
if isempty(mods)
    append_module(MODULES.mean_squared_error);   
end

% SVD removed 2013-10-24
%mods = find_modules(STACK, 'correlation', true);
%if isempty(mods)
%    append_module(MODULES.correlation);    
%end

META.perf_metric = @pm_nmse;
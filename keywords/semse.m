function semse()

global MODULES META STACK;

mods = find_modules(STACK, 'mean_squared_error', true);
if isempty(mods)
    append_module(MODULES.mean_squared_error.mdl(struct('norm_by_se', 1)));   
end

META.perf_metric = @pm_nmse;
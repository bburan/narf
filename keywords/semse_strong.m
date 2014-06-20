function semse_strong()

global MODULES META STACK;

mods = find_modules(STACK, 'mean_squared_error', true);
if isempty(mods)
    disp('semse_strong: 0.7');
    append_module(MODULES.mean_squared_error.mdl(struct('norm_by_se', 0.7)));   
end

META.perf_metric = @pm_nmse;

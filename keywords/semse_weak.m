% function semse_weak()
%
% weak shrinkage (norm_by_se=0.25)
function semse_weak()

global MODULES META STACK;

mods = find_modules(STACK, 'mean_squared_error', true);
if isempty(mods)
    disp('semse_weak: 0.4');
    append_module(MODULES.mean_squared_error.mdl(struct('norm_by_se', 0.4))); 
end

META.perf_metric = @pm_nmse;

function llpoisson()

global MODULES META STACK;

mods = find_modules(STACK, 'likelihood_poisson', true);
if isempty(mods)    
    append_module(MODULES.likelihood_poisson);    
end

META.perf_metric = @pm_nlogl;

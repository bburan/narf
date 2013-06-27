function [termcond, n_iters] = fit_anneal(options)

if nargin < 1
    options = saoptimset('MaxIter', 9000, ...
                         'MaxFunEvals', 9000, ...
                         'TolFun', 1e-12);  
end

[termcond, n_iters] = default_fitter_loop('simulannealbnd()', ...
    @(obj_fn, phi_init) simulannealbnd(obj_fn, phi_init, [], [], options));

function [term_cond, term_score, n_iters, term_step] = fit_anneal(options)

if nargin < 1
    options = saoptimset('MaxIter', 9000, ...
                         'MaxFunEvals', 9000, ...
                         'TolFun', 1e-12);  
end

[term_cond, term_score, n_iters, term_step] = default_fitter_loop('simulannealbnd()', ...
    @(obj_fn, phi_init) simulannealbnd(obj_fn, phi_init, [], [], options));

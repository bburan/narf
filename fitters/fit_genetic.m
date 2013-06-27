function [termcond, n_iters] = fit_genetic(options)

% TODO: help gaoptimset .... There are about a gazillion parameters to play with!
if nargin < 1
    options = gaoptimset('TolFun', 1e-12);  
end

[termcond, n_iters] = default_fitter_loop('ga()', ...
    @(obj_fn, phi_init) ga(obj_fn, length(phi_init), [], [], [], [], [], [], [], options));

end
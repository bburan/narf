function [termcond, n_iters, term_stepsize] = default_fitter_loop(fittername, highlevel_fn, bequiet)

global STACK META;

if ~exist('bequiet', 'var')
    bequiet = false;
end

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

n_iters = 1;
start_depth = find_fit_start_depth(STACK);
depths = find_fit_param_depths(STACK);
[mse, pen] = META.perf_metric();
score_prev = mse + pen;
prev_phi = [];

function score = my_obj_fn(phi)
    % Find the point at which the stack needs to be recalculated.
    if isempty(prev_phi)
        prev_phi = phi;
        unpack_fittables(phi);
        calc_xxx(start_depth);
    else
        idx_of_first_different_param = find(phi ~= prev_phi, 1);
        prev_phi = phi;
        unpack_fittables(phi);
        calc_xxx(depths(idx_of_first_different_param));
    end
    [m, p] = META.perf_metric();
    score = m + p;
    score_delta = score_prev - score;
    score_prev = score;
    
    % Print a newline and important info every 100 iterations
    if isequal(mod(n_iters, 100), 1) && ~bequiet
        fprintf('\n[Iter: %5d, score:%f delta:%6f]', n_iters, score, score_delta);
        % dbtickqueue(n_iters);
    end
    
    % Print 1 progress dot for every 5 iterations
    if isequal(mod(n_iters, 5), 0) %% && ~bequiet
        fprintf('.');
    end
    
    n_iters = n_iters + 1;
end

fprintf('----------------------------------------------------------------------\n');
fprintf('Fitting %d variables with %s\n', length(phi_init), fittername);

if nargout(highlevel_fn) ==4 
    [phi_best, ~, termcond, term_stepsize] = highlevel_fn(@my_obj_fn, phi_init);
else    
    [phi_best, ~, termcond] = highlevel_fn(@my_obj_fn, phi_init);
    term_stepsize = 1;
end

unpack_fittables(phi_best);

fprintf('Complete fit with %d objective function evaluations.\n', n_iters);
fprintf('----------------------------------------------------------------------\n');
end

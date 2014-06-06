function [term_cond, term_score, n_iters, term_step] = default_fitter_loop(fittername, highlevel_fn, bequiet)

global STACK META;

if ~exist('bequiet', 'var')
    bequiet = false;
end

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    term_cond = NaN;
    term_score = NaN;
    n_iters = 0;
    term_step = NaN;
    return 
end

n_iters = 1;
start_depth = find_fit_start_depth(STACK);
depths = find_fit_param_depths(STACK);
[mse, pen] = META.perf_metric();
score_prev = mse + pen;
prev_phi = [];

function [score, iters] = my_obj_fn(phi)
    iters = n_iters;
    
    % Find the point at which the stack needs to be recalculated.
    if isempty(prev_phi)
        % First time only
        prev_phi = phi;
        unpack_fittables(phi);
        update_xxx(start_depth);
    else
        % After that, use depths() to choose start point
        idx_of_first_different_param = find(phi ~= prev_phi, 1);
        prev_phi = phi;
        unpack_fittables(phi);
        update_xxx(depths(idx_of_first_different_param));
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
    
    % Print 1 progress dot for every 5 iterations no matter what
    if isequal(mod(n_iters, 5), 0) %% && ~bequiet
        fprintf('.');
    end
    
    n_iters = n_iters + 1;
end

fprintf('Fitting %d variables with %s\n', length(phi_init), fittername);

[term_phi, term_score, term_cond, term_step] = highlevel_fn(@my_obj_fn, phi_init);

unpack_fittables(term_phi);
idx_of_first_different_param = find(term_phi ~= prev_phi, 1);
update_xxx(depths(idx_of_first_different_param));

fprintf('Complete fit with %d objective function evaluations.\n', n_iters);
fprintf('----------------------------------------------------------------------\n');
end

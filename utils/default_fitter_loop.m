function [termcond, n_iters] = default_fitter_loop(fittername, highlevel_fn, bequiet)

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
[mse, pen] = META.perf_metric();
score_prev = mse + pen;

function score = my_obj_fn(phi)
    unpack_fittables(phi);
    calc_xxx(start_depth);
    [m, p] = META.perf_metric();
    score = m + p;
    score_delta = score_prev - score;
    score_prev = score;
    
    % Print a newline and important info every 100 iterations
    if isequal(mod(n_iters, 100), 1) && ~bequiet
        fprintf('\n[Iter: %5d, score:%f delta:%6f]', n_iters, score, score_delta);
    end
    
    % Print 1 progress dot for every 5 iterations
    if isequal(mod(n_iters, 5), 0) && ~bequiet
        fprintf('.');
    end
    
    n_iters = n_iters + 1;
end

fprintf('----------------------------------------------------------------------\n');
fprintf('Fitting %d variables with %s\n', length(phi_init), fittername);

[phi_best, ~, termcond] = highlevel_fn(@my_obj_fn, phi_init);
                                
unpack_fittables(phi_best);

fprintf('Complete fit in %d iterations.\n', n_iters);
fprintf('----------------------------------------------------------------------\n');
end

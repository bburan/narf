function [termcond] = default_fitter_loop(highlevel_fn, perf_metric)

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

if ~exist('iter_count', 'var')
    iter_count = 1;
end

start_depth = find_fit_start_depth(STACK);

function score = my_loop_fn(phi)
    % Print 1 progress dot for every 20 iterations, and newline every 1000
    if isequal(mod(cnt, 1000), 1)
        fprintf('\n[Iteration %d]', cnt);
    end
    if isequal(mod(cnt, 20), 1)
        fprintf('.');
    end
    cnt = cnt + 1;

    unpack_fittables(phi);
    calc_xxx(start_depth);
    score = XXX{end}.(objective_score);
end


fprintf('Fitting %d variables with fminunc()\n', length(phi_init));

len = length(flatten_field(XXX{end}.dat, XXX{end}.training_set, 'respavg'));

[phi_best, ~, termcond] = fminunc(@my_obj_fn, phi_init, options);
                                
unpack_fittables(phi_best);
end
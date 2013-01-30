function termcond = fit_fminunc(objective_score, options)

global XXX STACK;
cnt = 1;  % For printing progress dots

if nargin < 1
    objective_score = 'score';
end
if nargin < 2
    options = optimset('MaxIter', 9000, ...
                       'MaxFunEvals', 9000, ...
                       'TolFun', 1e-12, ...
                       'TolX', 1e-9);  
end

start_depth = find_fit_start_depth(STACK);

function score = my_obj_fn(phi)
    % Print 1 progress dot for every 20 iterations, and newline every 1000
    if isequal(mod(cnt, 1000), 1)
        fprintf('\n[Iteration %d]', cnt);
    end
    if isequal(mod(cnt, 20), 1)
        fprintf('.');
    end
    cnt = cnt + 1;

    unpack_fittables(phi);
    recalc_xxx(start_depth);
    score = XXX{end}.(objective_score);
end

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

recalc_xxx(1); 
fprintf('Fitting %d variables with fminunc()\n', length(phi_init));

len = length(flatten_field(XXX{end}.dat, XXX{1}.training_set, 'respavg'));

[phi_best, ~, termcond] = fminunc(@my_obj_fn, phi_init, options);
                                
unpack_fittables(phi_best);
end
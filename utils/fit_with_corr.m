function fit_with_corr()
% Fits any fittable fields using the training set correlation performance
% metric

global XXX STACK;
cnt = 0;  % For printing progress dots

    function thescore = my_scorer(phi, start_depth)
    % Perform the nonlinear fitting routine, starting at start_depth and
    % recalculating until the end of the XXX datastructure. 
    
    cnt = cnt + 1;
    % Print 1 progress dot for every 20 iterations, and newline every 1000
    if isequal(mod(cnt, 1000), 1)
        fprintf('\n[Iteration %d]', cnt);
    end
    if isequal(mod(cnt, 20), 1)
        fprintf('.');
    end
    
    unpack_fittables(phi);
    recalc_xxx(start_depth);
    thescore = XXX{end}.score_train_corr;
end

% The optimization
phi_init = pack_fittables(STACK);
fprintf('Fitting %d variables with fminsearch()\n', length(phi_init));

options = optimset('MaxIter', 1000, ...
                   'MaxFunEvals', 5000, ... % 100*numel(phi_init), ...
                   'TolFun', 1e-12, 'TolX', 1e-9);
phi_best = fminsearch(@my_scorer, phi_init, options);
unpack_fittables(phi_best);
end
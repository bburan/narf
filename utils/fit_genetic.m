function termcond = fit_genetic(objective_score, options)
% termcond = fit_genetic(objective_score, options)
%
% Fits all parameters in STACK marked with 'fit_fields' such that at the
% end of the STACK, the 'objective_score' field is minimized.
%
% ARGUMENTS:
%    objective_score    The name of the signal to use as objective score.
%                       Defaults to 'score' if no argument is passed

%
% RETURNS:
%    termcond    Termination condition of optimization.

global XXX STACK;
cnt = 1;  % For printing progress dots

if nargin < 1
    objective_score = 'score';
end
if nargin < 2
%     options = gaoptimset('MaxIter', 9000, ...
%                          'MaxFunEvals', 9000, ...
%                          'TolFun', 1e-12);  
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

fprintf('Fitting %d variables with genetic algorithm ga()\n', length(phi_init));

[phi_best, ~, termcond] = ga(@my_obj_fn, length(phi_init));
unpack_fittables(phi_best);

end
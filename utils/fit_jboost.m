function termcond = fit_jboost(objective_score)
% termcond = fit_jboost(objective_score, options)
% Boosts on 10 jackknifes, then shrinks
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

start_depth = find_fit_start_depth(STACK);

function score = my_obj_fn(phi)
    % Print 1 progress dot for every 20 iterations, and newline every 1000
%     if isequal(mod(cnt, 1000), 1)
%         fprintf('\n[Iteration %d]', cnt);
%     end
%     if isequal(mod(cnt, 20), 1)
%         fprintf('.');
%     end
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
fprintf('Fitting %d variables with Boosting()\n', length(phi_init));

[phi_best, score_best] = boosting(@my_obj_fn, phi_init', @(n,x,s) (n > 200), 1);
unpack_fittables(phi_best);
termcond = NaN;

end
function termcond = fit_boost(objective_score, iterations)
% termcond = fit_boost(objective_score, options)
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

if nargin < 1
    objective_score = 'score';
end
if nargin < 2
    iterations = length(pack_fittables(STACK)) * 10;
end

start_depth = find_fit_start_depth(STACK);

function score = my_obj_fn(phi)
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

fprintf('Fitting %d variables with fit_boost()\n', length(phi_init));

[phi_best, score_best] = boosting(@my_obj_fn, phi_init', @(n,x,s) (n > iterations), 1);
unpack_fittables(phi_best);
termcond = NaN;

end
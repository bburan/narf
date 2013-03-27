function n_steps_taken = fit_boost(objective_score, max_n_steps, min_stepsize, min_scoredelta)
% n_steps_taken = fit_boost(objective_score, n_steps, minstepsize, min_scoredelta)
%
% A generic NARF Fitting Routine which uses a boosting algorithm. Works
% well for sparse linear spaces initialized with magnitudes between
% negative one and positive one. Works more poorly outside that range, and
% will work very poorly for highly nonlinear systems.
%
% All arguments are optional. 
%
% ARGUMENTS:
%    objective_score    The name of the signal to use as objective score.
%                       Defaults to 'score' if no argument is passed
%
%    max_n_steps        Terminate search after this many steps are taken.
%                       Default: max(50, #parameters x 5)
%
%    min_stepsize       Terminate search if stepsize smaller than this. 
%                       Default: 10^-9
%
%    min_scoredelta     Terminate search if objective score improves less
%                       than this amount on each step.
%                       Default: 10^-12
%
% RETURNS:
%    n_steps            The number of boosting steps taken.

global XXX STACK;

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    n_steps_taken = 0;
    return 
end

if ~exist('objective_score','var'),  
    objective_score = 'score';
end

if ~exist('max_n_steps', 'var'),
    max_n_steps = max(50, length(phi_init(:)) * 5);
end

if ~exist('min_stepsize', 'var'),
    min_stepsize = 10^-9;
end

if ~exist('min_scoredelta', 'var'),
    min_scoredelta = 10^-12;
end

start_depth = find_fit_start_depth(STACK);

function score = my_obj_fn(phi)
    unpack_fittables(phi);
    recalc_xxx(start_depth);
    score = XXX{end}.(objective_score);
end

fprintf('Fitting %d variables with fit_boost()\n', length(phi_init));

[phi_best, ~, n_steps_taken] = boosting(@my_obj_fn, phi_init, ...
    @(n,x,s,d) (n > max_n_steps | s < min_stepsize | d < min_scoredelta), 1);

unpack_fittables(phi_best);

end
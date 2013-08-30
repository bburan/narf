function [term_cond, term_score,  n_iters, options] = fit_boo(varargin)
%  [termcond,n_iters, term_stepsize] = fit_boo(varargin)
%
% A clean interface to the boost algorithms. 
% 
%   'StopAtSeconds'        Stop after this many seconds.
%   'StopAtStepNumber'     Stop when the stepnumber is this small.
%   'StopAtStepSize'       Stop when stepsize is this small.
%   'StopAtRelScoreSize'   Stop when the relative improvement to the score
%                          is smaller than this level.
%   'StopAtAbsScoreSize'   Stop when the absolute improvement to the score 
%                          is smaller than this level.
%   'StepSize'             Initial step size.
%   'StepRel'              When true, step sizes are scaled by the variance 
%                          they create in the output.
%   'StepRelRecalcEvery'   Recalculate scaling factor every this many steps
%   'StepRelMin'           Minimum scaling factor bound for StepRel.
%   'StepRelMax'           Maximum scaling factor bound for StepRel.
%   'StepGrowth'           Multiply the next step size by this number if
%                          this step is successful. 
%   'StepShrink'           Multiply the next step size by this number if 
%                          no better step was found. 
% 
%   'Elitism'              When true, take a few steps searching only the
%                          most elite parameters (best effects on score).
%   'EliteParams'          Number of parameters to consider as elite.
%   'EliteSteps'           Number of steps to take on elite parameters
%                          before doing a wide search again.
% 
%  Finally, if you need truly arbitrary stopping conditions, you may set:
%    'TermFn'              When defined, this is used instead of the above
%                          termination conditions to allow for truly
%                          arbitrary stopping conditions.
%
% The termination function should accept 3 arguments:
%                     n    Step number of this iteration.
%                     s    Stepsize taken in this past step
%                     d    Score improvement (delta) vs previous point
%
% Boosting will terminate when TermFn returns a non-zero value.

global STACK;
phi = pack_fittables(STACK);
n_params = length(phi);

default_options = struct('StopAtSeconds', 60*60*24, ...
                         'StopAtStepNumber', max(50, n_params*5), ...
                         'StopAtStepSize', 10^-9, ...
                         'StopAtRelScoreDelta', 10^-12, ...
                         'StopAtAbsScoreDelta', 10^-12, ...
                         'InitStepSize', 1.0, ...
                         'StepRel', false, ...
                         'StepRelRecalcEvery', 1, ...
                         'StepRelMin', 10^-1, ...
                         'StepRelMax', 10^1, ...
                         'StepGrowth', 1.0, ...
                         'StepShrink', 0.5, ...
                         'Elitism', false, ...
                         'EliteParams', 10, ...
                         'EliteSteps', 1);

if (length(varargin) == 1)
    options = merge_structs(default_options, varargin{1});
else
    options = merge_structs(default_options, struct(varargin{:}));
end

if ~isfield(options, 'TermFn') 
    options.TermFn = create_term_fn(options);
end

[term_cond, term_score, n_iters, term_step] = ...
    default_fitter_loop('fit_boo()', ...
        @(obj_fn, phi_init) boost_algorithm(obj_fn, phi_init, options), true);

% If we are running this in an iterative fitter, next time start from here
options.InitStepSize = term_step;
end


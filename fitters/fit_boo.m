function [termcond, n_iters, term_stepsize] = fit_boo(options)
%  [termcond,n_iters, term_stepsize] = fit_boo(options)
%
% A consolidation of various boost algorithms that have been tried.
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
% The termination function which accepts 4 arguments:
%                     n    Step number of this iteration.
%                     x    The present x being considered.
%                     s    Stepsize taken in this past step
%                     d    Score improvement (delta) vs previous point
%
% Boosting will terminate when TermFn returns a true.

global STACK;
phi = pack_fittables(STACK);

n_params = length(phi);

default_options = struct('StopAtSeconds', 60*60*24, ...
                         'StopAtStepNumber', max(50, n_params*5), ...
                         'StopAtStepSize', 10^-9, ...
                         'StopAtRelScoreDelta', 10^-5, ...
                         'StopAtAbsScoreDelta', 10^-12, ...
                         'StepSize', 1.0, ...
                         'StepRel', false, ...
                         'StepRelRecalcEvery', 1, ...
                         'StepRelMin', 10^-3, ...
                         'StepRelMax', 10^3, ...
                         'StepGrowth', 1.0, ...
                         'StepShrink', 0.5, ...
                         'Elitism', false, ...
                         'EliteParams', 10, ...
                         'EliteSteps', 1);

if ~exist('options', 'var')
    opts = default_options;
else
    opts = merge_structs(default_options,options);
end

if ~opts.Elitism
    opts.EliteParams = n_params;
    opts.EliteSteps = 1;
end

t_start = clock(); %
first_5_deltas = [];       % Used for relative stopping

function stopcode = termfn(n,x,s,d)
    if (n > opts.StopAtStepNumber),   
        stopcode = 1; 
        fprintf('Took more than %d steps in total, stopping.\n', ...
            opts.StopAtStepNumber);
        return;
    end
    if (s < opts.StopAtStepSize),      
        stopcode = 2;  
        fprintf('Stepsize less than absolute limit %f, stopping.\n', ...
            opts.StopAtStepSize);
        return; 
    end
    if (d < opts.StopAtAbsScoreDelta), 
        stopcode = 3;
        fprintf('Absolute score improvement less %f, stopping.\n', ...
            opts.StopAtAbsScoreDelta);
        return; 
    end
    
    if length(first_5_deltas) < 5
        first_5_deltas(end+1) = d;
    else
        if (d / mean(first_5_deltas(2:end))) < opts.StopAtRelScoreDelta
            stopcode = 4;
            fprintf('Relative score improvement less %f, stopping.\n', ...
                opts.StopAtRelScoreDelta);
            return
        end
    end
       
    if (etime(clock, t_start) > opts.StopAtSeconds)
        stopcode = 5; 
        fprintf('Fit time exceeded %f, stopping.\n', ...
                opts.StopAtSeconds);
        return
    end

    stopcode = false;   
end

if ~isfield(opts, 'TermFn')
    opts.TermFn = @termfn;
end

[termcond, n_iters, term_stepsize] = default_fitter_loop('fit_boo()', ...
    @(obj_fn, phi_init) boost_algorithm(obj_fn, phi_init, opts), true);

end


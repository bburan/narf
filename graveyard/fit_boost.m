function [termcond, n_iters] = fit_boost(max_n_steps, min_stepsize, min_scoredelta, relative_delta, vary_stepsize, starting_stepsize)
% [termcond, n_iters] = fit_boost(n_steps, minstepsize, min_scoredelta, relative_delta, vary_stepsize, starting_stepsize)
%
% A generic NARF Fitting Routine which uses a boosting algorithm. Works
% well for sparse linear spaces initialized with magnitudes between
% negative one and positive one. Works more poorly outside that range, and
% will work very poorly for highly nonlinear systems.
%
% All arguments are optional. 
%
% ARGUMENTS:
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
%    relative_delta     When true, the min_scoredelta will act relative to
%                       the average of the first five delta steps. 
%
%    vary_stepsize      When true, stepsizes will vary to make them all
%                       have roughly the same effect on the score
%
% RETURNS:
%    n_steps            The number of boosting steps taken.

global STACK;

phi_init = pack_fittables(STACK);

if ~exist('max_n_steps', 'var'),
    max_n_steps = max(50, length(phi_init(:)) * 5);
end

if ~exist('min_stepsize', 'var'),
    min_stepsize = 10^-9;
end

if ~exist('min_scoredelta', 'var'),
    min_scoredelta = 10^-12;
end

if ~exist('relative_delta', 'var'),
    relative_delta = false;
end

if ~exist('vary_stepsize', 'var'),
    vary_stepsize = false;
end
if ~exist('starting_stepsize', 'var'),
    starting_stepsize = 1;
end



first5_deltas = nan(5,1); % First delta never used

function stop = term_fn(n,x,s,d)    
    if (n > max_n_steps)
        stop = 1
        return
    end
    if (s < min_stepsize)
        stop = 2
        return
    end
    if (d < min_scoredelta) && ~relative_delta
        stop = 3
        return
    end
    if relative_delta        
        idx = find(isnan(first5_deltas), 1);
        if ~isempty(idx)
            first5_deltas(idx) = d;
        else
           if (d / mean(first5_deltas(2:end))) < min_scoredelta
                stop = 4
                return
           end
        end
    end
    stop = false;   
end

[termcond, n_iters] = default_fitter_loop('fit_boost()', ...
    @(obj_fn, phi_init) boosting(obj_fn, phi_init, @term_fn, starting_stepsize, 2, vary_stepsize), true);    

end
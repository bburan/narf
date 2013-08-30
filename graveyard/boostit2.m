function [termcond, n_iters] = boostit2()
% Boosting algorithm which:
% 1. Fits modules one at a time.
% 2. Boosts N steps, where N = # of fittable parameters in that module.
%    The initial step size is equal to the size of the largest free param.
% 3. Iterates, fitting every module 20*M times, where M is the # of modules
% 4. Does NOT early stop.

global STACK;

N_iters = 10 * sum(cellfun(@(a) sum(cellfun(@(b) isfield(b, 'fit_fields') && ~isempty(b.fit_fields), a)), STACK));

min_scoredelta = 10^-4;
max_n_steps = 1000000; % (Will be overwritten but must start large)
steps_taken = 0;
min_stepsize = 10^-4;

stepsize_cache = [];

function [termcond, n_iters] = fitter_fn()
    % Yes, this looks goofy, but it has to be here for scope reasons
    % and so that steps_taken is updated properly.
    function stop = term_fn(n,x,s,d)
        steps_taken = n;       

        if (n > max_n_steps)
            stop = 1;
            return;
        end
    
        if (s < min_stepsize)
            stop = 2;
            return;
        end
        
        if (d < min_scoredelta)
            stop = 3;
            return;
        end        

        stop = false;   
    end
    
    phi = pack_fittables(STACK);
    initial_step_size = max(phi(:));
    if (initial_step_size <= 0)
        initial_step_size = 1;
    end
    max_n_steps = length(phi(:)) + steps_taken;
    [termcond, n_iters] = default_fitter_loop('fit_boost()', ...
        @(obj_fn, phi_init) boosting(obj_fn, phi_init, @term_fn, initial_step_size), true);
end

[termcond, n_iters] = fit_iteratively(@fitter_fn, N_iters);

end
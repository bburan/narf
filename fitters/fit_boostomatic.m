function [termcond, n_iters] = fit_boostomatic(max_n_steps, min_stepsize, min_scoredelta, growthfactor, num_elites, num_elite_iters)

global STACK;

phi_init = pack_fittables(STACK);

if ~exist('max_n_steps', 'var'),
    max_n_steps = max(50, length(phi_init(:)) * 10);
end

if ~exist('min_stepsize', 'var'),
    min_stepsize = 10^-6;
end

if ~exist('min_scoredelta', 'var'),
    min_scoredelta = 10^-12;
end

if ~exist('growthfactor', 'var'),
    growthfactor = sqrt(2);
end

if ~exist('num_elites', 'var'),
    num_elites = 1;
end

if ~exist('num_elite_iters', 'var'),
    num_elite_iters = 5;
end


function stop = term_fn(n,x,s,d)    
    fprintf('Step [%d/%d]\n', n, max_n_steps);
    if (n > max_n_steps)
        stop = 1;
        return
    end
    %if (s < min_stepsize)
    %    stop = 2;
    %    return
    %end
    %if (d < min_scoredelta)
    %    stop = 3;
    %    return
    %end

    stop = false;   
end

[termcond, n_iters] = default_fitter_loop('fit_boostomatic()', ...
    @(obj_fn, phi_init) boostomatic(obj_fn, phi_init, @term_fn, growthfactor, num_elites, num_elite_iters), true);    

end
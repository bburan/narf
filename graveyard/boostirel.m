function boostirel()
% Iterative boosting with globalized early stopping
% Stops each iteration when relative improvement to GLOBAL performance is
% reaching 0.1% or less of the average of the first 4-5 steps

earlystop = 10^-3; % Try 10^-3 to 10^-7 or so?

global STACK;
phi_init = pack_fittables(STACK);
max_n_steps = max(50, length(phi_init(:)) * 5);
min_stepsize = 10^-9;
min_scoredelta = earlystop;
first5_deltas = nan(5,1); % First delta never used

function stop = term_fn(n,x,s,d)    
    if (n > max_n_steps) || (s < min_stepsize)
        stop = 1;
        return;
    end
        
	idx = find(isnan(first5_deltas), 1);
    if ~isempty(idx)
        first5_deltas(idx) = d;
    else
       if (d / mean(first5_deltas(2:end))) < min_scoredelta
            stop = 4;
            return
       end
    end     
    stop = false;   
end

fitter = @() default_fitter_loop('fit_boost()', ...
    @(obj_fn, phi_init) boosting(obj_fn, phi_init, @term_fn, 1), true);    

fit_iteratively(fitter, 10);

end
function mixfit3()

global STACK;

mse();


% alternative to qboost, run with quicker stop limit
phi_init = pack_fittables(STACK);
%fit_boost(n_steps, minstepsize, min_scoredelta, relative_delta)
fit_boost(length(phi_init),10^-9,10^-2,true,true);

% iterative part (boostirelsvd() modified for variable delta)

earlystop = 10^-4; % Try 10^-3 to 10^-7 or so?

phi_init = pack_fittables(STACK);
max_total_steps=max(50, length(phi_init(:)) * 2)*2;
max_n_steps = 10;
max_loops=ceil(max_total_steps./max_n_steps);
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
    @(obj_fn, phi_init) boosting(obj_fn, phi_init, @term_fn, 1, 10, ...
                                 true), true);    

fit_iteratively(fitter, max_loops);


% final qboost with variable stepsize
fit_boost(length(phi_init),10^-9,10^-3,true,true);

end
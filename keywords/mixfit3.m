function mixfit3(SKIP_MSE)

global STACK;

if ~exist('SKIP_MSE','var'),
    SKIP_MSE=0;
end
if ~SKIP_MSE,
    mse();
end

% alternative to qboost, run with quicker stop limit
min_stepsize = 10^-7;
min_scoredelta = 10^-5; % Try 10^-3 to 10^-7 or so?;
phi_init = pack_fittables(STACK);

%fit_boost(n_steps, minstepsize, min_scoredelta, relative_delta, vary_stepsize, starting_stepsize)
fit_boost(length(phi_init),min_stepsize,10^-2,true,true);

% iterative part (boostirelsvd() modified for variable delta)
phi_len_max=0;
for ii=1:length(STACK),
    mm = STACK{ii}{1};
    if isfield(mm,'fit_fields'),
        phi_len=0;
        for jj = 1:length(mm.fit_fields),
            p = mm.fit_fields{jj};
            phi_len=phi_len+numel(mm.(p));
        end
        phi_len_max=max(phi_len_max,phi_len);
    end
end
max_total_steps=max(25, phi_len_max) * 4;
max_n_steps = 10;
max_loops=ceil(max_total_steps./max_n_steps);

first5_deltas = nan(5,1); % First delta never used
phi_init = pack_fittables(STACK);

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
fit_boost(length(phi_init),min_stepsize,min_scoredelta,true,true);

end
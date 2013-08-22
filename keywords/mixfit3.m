function mixfit3(SKIP_MSE)

global STACK;

if ~exist('SKIP_MSE','var'),
    SKIP_MSE=0;
end
if ~SKIP_MSE,
    mse();
end

min_stepsize = 10^-7; % stop boosting if abs step size smaller than this
min_scoredelta = 10^-5; % stop boosting if mse improvement smaller than this

%% initially run with qboost with large stop limit
phi_init = pack_fittables(STACK);
%fit_boost(n_steps, minstepsize, min_scoredelta, relative_delta, vary_stepsize, starting_stepsize)
fit_boost(length(phi_init),min_stepsize,10^-2,true,true);


%% iterative part

% figure out max number of parameters
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

% use phi_len_max to figure out how many iterations to complete
max_steps_per_iteration=min(10,phi_len_max);
iteration_frac=max_steps_per_iteration./phi_len_max;
max_loops=round(10./iteration_frac);

first5_deltas = nan(5,1); % First delta never used

function stop = term_fn(n,x,s,d)
    
    % iterations per loop depends on the number of parameters to fit
    max_n_steps=min(10,length(x));
    
    if (n >= max_n_steps)
        fprintf('%d steps completed, stopping\n',max_n_steps);
        stop = 1;
        return;
    end
    if (s < min_stepsize)
        fprintf('below min step size, stopping\n');
        stop = 2;
        return;
    end
     if d<min_scoredelta,
        fprintf('score delta %.8f less than %.8f, stopping\n',...
                d,min_scoredelta);
        stop=3;
        return
    end
    stop = false;
end

phi_init = pack_fittables(STACK);
fitter = @() default_fitter_loop('fit_boost()', ...
    @(obj_fn, phi_init) boosting(obj_fn, phi_init, @term_fn, 1, 2, ...
                                 true), true);

fit_iteratively(fitter, max_loops);

% final qboost with variable stepsize
fit_boost(length(phi_init),min_stepsize,min_scoredelta,true,true);

end
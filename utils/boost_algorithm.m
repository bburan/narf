function [x_bst, s_bst, termcond] = boost_algorithm(objfn, x_0, options)
%  [x_bst, s_bst, termcond] = boost_algorithm(objfn, x_0, options)
%
% See documentation of fit_boo.m for more information.
%
global XXX;

n_params = length(x_0);

n = 0;             % Step number
x = x_0(:);        % Current search location
s = objfn(x);      % Current score
s_delta = s;       % Improvement of score over the previous step
stepsize = options.StepSize;  % Starting step size.
stim = flatten_field(XXX{end}.dat, XXX{end}.training_set, 'stim');       
    
effect    = zeros(n_params, 1);  % The variance it creates on the STIM 
eliteness = zeros(n_params, 1);  % The effect on the score of each parameter

last_elite_recalc = -100;
last_effect_recalc = -2*options.StepRelRecalcEvery;

while (true) 
    x_next = x;
    s_next = s;  
        
    if (options.Elitism && last_elite_recalc ~= n) || ...
       (options.StepRel && last_effect_recalc <= n - options.StepRelRecalcEvery),
   
        fprintf('Calculating parameter importance');        
        for d = 1:n_params
            stepdir = zeros(n_params, 1);
            stepdir(d) = stepsize;            
            
            % this next line actually does the work, since it causes a
            % calc_xxx to occur 
            eliteness(d) = objfn(x + stepdir) - s;
            
            if options.StepRel
                newstim = flatten_field(XXX{end}.dat, XXX{end}.training_set, 'stim');
                effect(d) = sum((stim - newstim).^2);    
            end            
        end
        
        if options.StepRel
            % Calculate effects such that the step taken will satisfy
            % the constraint: StepRelMin < step <  StepRelMax            
            effect = (effect./max(effect));  %% Normalize from 0 to 1                        
            cutoff = (options.StepRelMin / options.StepRelMax);   
            effect(effect <= cutoff) = cutoff; 
            effect = effect ./ options.StepRelMax;                    
            
            last_effect_recalc = n;
        end

        if options.Elitism
            idxs = 1:n_params;
            D = [eliteness, idxs'];
            D = sortrows(D);
            elites = D(1:options.EliteParams, 2);
            %fprintf('\nElite indexes: ');
            %fprintf('%d ', elites);
            %fprintf('\n');
            last_elite_recalc = n;
        end
        fprintf('\n');    
    end

    for iter = 1:options.EliteSteps
        fprintf('Boosting'); 
        for eidx = 1:options.EliteParams, 
            stepdir = zeros(n_params, 1);
            if options.Elitism
                idx = elites(eidx);
            else
                idx = eidx;
            end
            
            if options.StepRel
                stepdir(idx) = stepsize / effect(idx);
            else
                stepdir(idx) = stepsize;
            end 

            x_fwd = x + stepdir;
            x_bck = x - stepdir;            
            
            s_fwd = objfn(x_fwd);
            s_bck = objfn(x_bck);           
            
            % Prepare to step forward if that is better
            if s_fwd < s_next
                s_next = s_fwd;
                x_next = x_fwd;
            end
            
            % Prepare to step backward if that is better
            if s_bck < s_next
                s_next = s_bck;
                x_next = x_bck;
            end
            
        end    
                        
        if all(x == x_next)
            % Stepsize was too big so no step was taken
            stepsize = stepsize * options.StepShrink;
            fprintf('Stepsize shrunk to %d\n', stepsize);            
        else
            % Step was taken successful
            s_delta = s - s_next;   % Improvement in score            
            dirs = 1:n_params;      % Direction of the step;
            dir = dirs(x ~= x_next);
            
            x = x_next;   % Take a step in the best direction
            s = objfn(x); % Recalculate the score
            stim = flatten_field(XXX{end}.dat, XXX{end}.training_set, 'stim');       
    
            % Print the improvement after stepping
            fprintf('step: %d, coef# %3d, delta: %d, score:%d\n', stepsize, dir, s_delta, s);

            % Prepare for next loop
            stepsize = stepsize * options.StepGrowth;
            n = n + 1;
        end

        % Should we stop searching?
        if options.TermFn(n, x, stepsize, s_delta)
            x_bst = x;
            s_bst = s;
            termcond = options.TermFn(n, x, stepsize, s_delta);
            return;
        end
        
        % Force a recalc because we didn't take a step and are elite
        if all(x == x_next) && (options.Elitism)
           last_elite_recalc = n-100;  % Force another recalc
           break;
        end
    end        
end

end %% end function


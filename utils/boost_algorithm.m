function [term_phi, term_score, term_cond, term_step] = boost_algorithm(objfn, x_0, options)
%  [term_phi, term_score, term_cond, term_step] = boost_algorithm(objfn, x_0, options)
%
% See documentation of fit_boo.m for more information.
%
global XXX STACK GA_LowerBound GA_UpperBound GA_Bounded; % NARF_DEBUG NARF_DEBUG_FIGURE;

n_params = length(x_0);

if ~options.Elitism
    options.EliteParams = n_params;
    options.EliteSteps = 1;
end

n = 1;             % Step number
x = x_0(:);        % Current search location
s = objfn(x);      % Current score
s_delta = inf;     % Improvement of score over the previous step
stepsize = options.InitStepSize;  % Starting step size.

effect    = zeros(n_params, 1);  % The variance it creates on the STIM 
eliteness = zeros(n_params, 1);  % The effect on the score of each parameter

last_elite_recalc = -100;
last_effect_recalc = -2*options.StepRelRecalcEvery;

fprintf('\n');

% the following block handles the definition of the step scaling factors
steps_scaling_factor = ones(n_params,1);
if options.ScaleStepSize,
    initialize_GA_GLOBALS(x_0);
    for ii = 1:length(GA_Bounded)
        if GA_Bounded(ii),
            steps_scaling_factor(ii) = GA_UpperBound(ii) - GA_LowerBound(ii);
        end
    end
    % we scale the bounds so that their mediaan is equal to 1 (as to make
    % the inclusion of bounds as independent as possible from the stopping
    % step size)
    steps_scaling_factor = steps_scaling_factor / median(steps_scaling_factor);
end

while (true) 
    x_next = x;
    s_next = s;  
        
    if (options.Elitism && last_elite_recalc ~= n) || ...
       (options.StepRel && last_effect_recalc <= n - options.StepRelRecalcEvery),
           
        [~, o] = objfn(x); % Calculate base stim rate
        stim = flatten_field(XXX{end}.dat, XXX{end}.training_set, 'stim');       
    
        fprintf('Calculating effect deltas for small epsilon step');        
        for d = 1:n_params
            stepdir = zeros(n_params, 1);
            stepdir(d) = stepsize ./ 1000 ./ steps_scaling_factor(d); %unused?
            
            [s_del, o] = objfn(x);
            eliteness(d) = abs(s - s_del);
            
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
                stepdir(idx) = stepsize / effect(idx) / steps_scaling_factor(idx);
            else
                stepdir(idx) = stepsize / steps_scaling_factor(idx);
            end 

            x_fwd = x + stepdir;
            x_bck = x - stepdir;            
            
            [s_fwd, o] = objfn(x_fwd);
            [s_bck, o] = objfn(x_bck);           
            
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
            fprintf('stepsize shrunk to %d\n', stepsize);      
 
            % Should we stop searching? (perhaps stepsize is too small)
            term_cond = options.TermFn(n, stepsize, s_delta, o);
            if term_cond
                [s, o] = objfn(x);  % don't leave us in wrong state from last step
                term_phi = x;
                term_score = s;
                term_step = stepsize;
                return;
            end     
        else
            % Step was taken successful
            s_delta = s - s_next;   % Improvement in score            
            dirs = 1:n_params;      % Direction of the step;
            dir = dirs(x ~= x_next);
            
            actual_stepsize = sum(x_next-x);
            
            % Should we stop searching? (perhaps s_delta is too small)
            term_cond = options.TermFn(n, stepsize, s_delta, o);
            if term_cond
                [s, o] = objfn(x); 
                term_phi = x;
                term_score = s;
                term_step = stepsize;
                return;
            end
                       
            x = x_next; % Take a step in the best direction
            s = s_next; 
            
            % Print the improvement after stepping
            fprintf('eval: %d, size: %.3e, act: %.3e, coef#%3d, s_delta: %.3e, score:%e\n', ...
                     o, stepsize, actual_stepsize, dir, s_delta, s);

%             % Possibly display debug info as well
%             if exist('NARF_DEBUG', 'var') && NARF_DEBUG
%                 if isempty(NARF_DEBUG_FIGURE),
%                     NARF_DEBUG_FIGURE=figure;
%                 else
%                     sfigure(NARF_DEBUG_FIGURE);
%                 end
%                 if n==0,
%                     x0=x;
%                 end
%                 clf;
%                 plot(1:length(x), x0, 'k+--');
%                 hold on
%                 errorbar(1:length(x), x, stepsize./effect);
%                 plot(dir, x_next(dir), 'ro');
%                 hold off
%                 xlabel('parameter');
%                 drawnow
%             end
            
            % Prepare for next loop
            stepsize = stepsize * options.StepGrowth;
            n = n + 1;
        end
        
        % Force a recalc because we didn't take a step and are elite
        if all(x == x_next) && (options.Elitism)
           last_elite_recalc = n-100;  % Force another recalc
           break;
        end
    end        
end


end %% end function


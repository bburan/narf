function [x_bst, s_bst, termcond] = boostomatic(objfn, x_0, termfn, ...
        growth_rate, num_elite_params, elite_iterations)
% Basically like "boosting.m", but tries to be smarter
%     - Stepsizes can grow arbitrarily large since they grow continuously
%     - More time is spent boosting "useful" parameters
%     - All parameters are treated equally and may become elite at any time
%
% Num_elite_params: How many parameters will be boosted as 'elite'
% elite_iterations: How many times we should just boost on those elite
% parameters before going back to the global search.

if nargin ~= 6,  
    error('boostomatic() needs 6 arguments.');
end

stepsize = 1;  % Starts at 1 but will grow
stepscale = 2; % Halve the size if there is a problem

% Good defaults?
% growth_rate = 1.1;

% Starting search point
x = x_0(:);
s = objfn(x);
l = length(x_0);
n = 0;  % Step number
s_delta = s;

eliteness = zeros(size(x_0));

while (true) 
    x_next = x;  % x_next holds best stepping direction so far
    s_next = s;
    
    % We cannot assume that the performance metric will be MSE in every
    % case; therefore you should query META.perf_metric() directly or pass
    % it from some function that is already calling perf_metric(). We'll
    % take the former case to avoid complicating the interface.    
    global META;
    [~,~,val_s] = META.perf_metric();
    
    fprintf('pm_est:  %12.8f  pm_val:  %12.8f\n', s, val_s);     
    fprintf('Wide Search');
    
    % Estimate the most elite parameters 
    for d = 1:l
        stepdir = zeros(l, 1);
        stepdir(d) = stepsize;
                
        % Step in the direction
        x_fwd = x + stepdir;
        x_bck = x - stepdir;

        s_fwd = objfn(x_fwd);
        s_bck = objfn(x_bck);
        
        eliteness(d) = min(s_fwd, s_bck); 
        fprintf('.');
    end
    
    idxs = 1:length(x_0);
    D = [eliteness(:), idxs'];
    D = sortrows(D);
    elites = D(1:num_elite_params, 2);
    fprintf('\nElite indexes:');
    fprintf('%d ', elites);
    fprintf('\n');
    
    for iter = 1:elite_iterations               
        fprintf('Elitist Search');    
        for d = 1:min(l, num_elite_params)    
            stepdir = zeros(l, 1);
            stepdir(elites(d)) = stepsize;
            
            % Step in the direction
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
            
            % Print a 20 dot progress bar
            if mod(d, ceil(l/20)) == 1
                fprintf('.');
            end
        end    
                        
        % If the search point has not changed, step size is too big.
        if all(x == x_next)
            stepsize = stepsize / stepscale;
            fprintf('Stepsize shrunk to %d\n', stepsize);
        else
            % Grow stepsize continuously
            stepsize = stepsize * growth_rate;
            
            % Compute the improvement in score
            s_delta = s - s_next;
            
            % Compute the direction of the step;
            dirs = 1:l;
            dir = dirs(x ~= x_next);
            
            % Take a step in the best direction
            x = x_next;
            % svd -- recalc error to reflect x_next in XXX
            s = objfn(x);
            
            % Print the improvement after stepping
            fprintf('coef# %3d, delta: %d, score:%d\n', dir, s_delta, s);
            n = n + 1;
            
            % Should we quit?
            if termfn(n, x, stepsize, s_delta)                
                % Return the best value found so far
                x_bst = x;
                s_bst = s;
                termcond = termfn(n, x, stepsize, s_delta);                
                return;                
            end
        end
            
    end
        
end

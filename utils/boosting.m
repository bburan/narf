function [x_bst, s_bst, termcond] = boosting(objfn, x_0, termfn, stepsize, stepscale, vary_stepsize)
% [x_bst, s_bst, termcond] = boosting(objfn, x_0, termfn, stepsize=1, stepscale=10, vary_stepsize=false)
%
% A naive 'boosting' search method in which stepsize can only decrease from
% its initial value. Good for linear searches of FIR coefficient space, but
% not suitable for nonlinear searches in which. 
%
% ARGUMENTS: 
%    objfn     Objective function. Must accept vector x and return a scalar
%              which is to be minimized.
%
%    x_0       Starting point for the optimization of vector x.
%
%    termfn    Termination function which accepts 4 arguments:
%                     n    Step number of this iteration.
%                     x    The present x being considered.
%                     s    Stepsize taken in this past step
%                     d    Score improvement (delta) vs previous point
%              The loop terminates when this returns a true.                
%
%    stepsize  Initial size of steps to take when optimizing x.
%              Default is 1. 
%
%    stepscale When no better steps are found, stepsize is scaled by:
%                stepsize = stepsize / stepscale.
%              Default is stepscale=2, the classic "binary search"
% 
% RETURNS:
%    x_bst     The best vector x found so far.
%
%    s_bst     The score of that vector, as evaluated by objfn.
%
%    termcond  The termination condition 
%
% DETAILS:
% Given a starting vector x = x_0, this algorithm samples a step in each
% direction along every dimension in x and moves to that new point 
% if it improves the score according to objective function objfn(x) more
% than any alternative step along a different dimension of x.
%
% Will terminate whenever @termfn(n,x,s,d) returns a boolean true value.
% Good termination function examples include:
%    @(n,x,s,d) n>300          Terminate after 300 steps.
%    @(n,x,s,d) d<10^-6        Terminate when is score delta is improving
%                              very slowly.
%    @(n,x,s,d) s<10^-3        Terminate when step size is small.

if nargin < 3,  
    error('boosting() needs 3 or more arguments.');
end

if ~exist('stepsize','var'),    
    stepsize = 1;
end

if ~exist('stepscale','var'),    
    stepscale = 2;
end
    
if (stepsize <= 0) 
    error('stepsize must be > 0');
end

if (stepscale <= 1) 
    error('stepscalemust be > 1');
end

if ~exist('vary_stepsize','var')
    vary_stepsize = false;   
end

global XXX;
global META;

% Starting search point
x = x_0(:);
s = objfn(x);
l = length(x_0);
n = 0;  % Step number
s_delta = s;
lastdeltastep=-1;
while ~termfn(n, x, stepsize, s_delta)    
    x_next = x;  % x_next holds best stepping direction so far
    s_next = s;
    
    deltas = ones(l, 1);
    stim = flatten_field(XXX{end}.dat, XXX{end}.training_set, 'stim');
    
    % if vary_stepsize, update delta every 5 steps
    if (vary_stepsize && ~mod(n,5) && n>lastdeltastep)  
        fprintf('Calculating Deltas for small Epsilon (stepsize/1000)');
        for d = 1:l
            stepdir = zeros(l, 1);
            % SVD alternative. epsilon=stepsize./1000;
            stepdir(d) = stepsize./1000;
            %stepdir(d) = x(d) * stepsize;             
            deltas(d) = abs(s - objfn(x + stepdir)); % To ensure evaluation
            newstim = flatten_field(XXX{end}.dat, XXX{end}.training_set,'stim');
            deltas(d) = sum((stim - newstim).^2);
    
            % Print a 20 dot progress bar
            if mod(d, ceil(l/20)) == 1
                fprintf('.');
            end
        end
        
        % SVD alternative, scale max delta to 1000 and truncate min
        % to 1.
        deltas=deltas./max(deltas).*10^3;
        deltas(deltas <= 1) = 1;
        
        % renormalize just in case spread of effects does not span 10^3
        deltas=deltas./min(deltas)./10; 
        
        % If there was a very minor change in the output, set it to 1 so
        % you don't see a huge explosion later.
        %deltas(deltas <= 10^-6) = 1;
        minidx=find(deltas==min(deltas),1);
        maxidx=find(deltas==max(deltas),1);
        
        fprintf('delta min %.3f (%d) max %.3f (%d)\n',...
                deltas(minidx),minidx,deltas(maxidx),maxidx);
        
        % record so we don't recalc deltas over and over for
        % stepsize reductions
        lastdeltastep=n;
        
        dbtickqueue(s);
    end
    
    % Try to take a step along every dimension
    fprintf('Boosting');
    
    for d = 1:l
        stepdir = zeros(l, 1);
        %stepdir(d) = stepsize + (stepsize/deltas(d));
        stepdir(d) = stepsize/deltas(d);
        if stepdir(d) == 0
            stepdir(d) = stepsize;
        end

        % Step in the direction
        x_fwd = x + stepdir;
        x_bck = x - stepdir;

        s_fwd = objfn(x_fwd);
        s_bck = objfn(x_bck);
        
        % Possibly record a step forward as being better
        if s_fwd < s_next
            s_next = s_fwd;
            x_next = x_fwd;
        end

        % Possibly record a backward step as being better
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
        fprintf('Decreased stepsize to %d\n', stepsize);
    else        
        % Compute the improvement in score
        s_delta = s - s_next;
        
        % Compute the direction of the step;
        dirs = 1:l;
        dir = dirs(x ~= x_next);
        
        % Take a step in the best direction
        x = x_next; 
        s = objfn(x);
        
        % Print the improvement after stepping
        fprintf('coef# %3d, delta: %d, score:%d\n', dir, s_delta, s);
        n = n + 1;
        
     end
     
     % We cannot assume that the performance metric will be MSE in every
     % case; therefore you should query META.perf_metric() directly or pass
     % it from some function that is already calling perf_metric(). We'll
     % take the former case to avoid complicating the interface.    
     if ~isempty(XXX{end}.test_set),
         [~,~,val_s] = META.perf_metric();
         fprintf('pm_est:  %12.8f  pm_val:  %12.8f\n', s, val_s);
     else
         fprintf('pm_est:  %12.8f\n', s);
     end
end

% Return the best value found so far
x_bst = x;
s_bst = s;
termcond = termfn(n, x, stepsize, s_delta);

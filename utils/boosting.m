function [x_bst, s_bst, n] = boosting(objfn, x_0, termfn, stepsize, stepscale)
% [x_bst, s_bst, n] = boosting(objfn, x_0, termfn, stepsize=1, stepscale=10)
%
% A naive 'boosting' search method in which stepsize can only decrease from
% its initial value. Good for linear searches of FIR coefficient space, but
% not suitable for nonlinear searches in which. 
%
% ARGUMENTS: 
%    objfn     Objective function. Must accept vector x and return a scalar
%              which is to be minimized.
%    x_0       Starting point for the optimization of vector x.
%    termfn    Termination function which accepts 4 arguments:
%                     n    Step number of this iteration.
%                     x    The present x being considered.
%                     s    Stepsize taken in this past step
%                     d    Score improvement (delta) vs previous point
%              The loop terminates when this returns a true.                
%    stepsize  Initial size of steps to take when optimizing x.
%    stepscale When no better steps are found, stepsize is scaled by:
%                stepsize = stepsize / stepscale.
%              Setting stepscale=2 results in classic "binary search"
% 
% RETURNS:
%    x_bst     The best vector x found so far.
%    s_bst     The score of that vector, as evaluated by objfn.
%    n         The number of boosting steps taken.
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
    stepscale = 10;
end
    
if (stepsize <= 0) 
    error('stepsize must be > 0');
end

if (stepscale <= 1) 
    error('stepscalemust be > 1');
end

% Starting search point
x = x_0(:);
s = objfn(x);
l = length(x_0);
n = 0;  % Step number
s_delta = s;

while ~termfn(n, x, stepsize, s_delta)    
    x_next = x;  % x_next holds best stepping direction so far
    s_next = s;
    
    % Try to take a step along every dimension
    fprintf('boosting.m: Stepping');
    
    for d = 1:l
        stepdir = zeros(l, 1);
        stepdir(d) = stepsize;
                
        % Step in the direction
        x_fwd = x + stepdir;
        x_bck = x - stepdir;

        s_fwd = objfn(x_fwd);
        s_bck = objfn(x_bck);
        
        % Take a step forward if that is better
        if s_fwd < s_next
            s_next = s_fwd;
            x_next = x_fwd;
        end

        % Take a step backward if that is better
        if s_bck < s_next
            s_next = s_bck;
            x_next = x_bck;
        end
        
        % Print a 20 dot progress bar
        if mod(d, ceil(l/20)) == 1
            fprintf('.');
        end
    end
    
    % Take a step in the best direction
    x = x_next; 
    s = s_next;
    
    % If the search point has not changed, step size is too big.
     if all(x == x_next)
         stepsize = stepsize / stepscale;
         fprintf('Stepsize decreased to %d\n', stepsize);
     else        
        % Compute the improvement in score
        s_delta = s - s_next;
    
        % Compute the direction of the step;
        dirs = 1:l;
        dir = dirs(x ~= x_next);
        
        
        % Take a step in the best direction
        x = x_next; 
        s = s_next;
        
        % Print the improvement after stepping
        fprintf('Step on %d improved by: %d\n', dir, s_delta);
        n = n + 1;
     end
    
end

% Return the best value found so far
x_bst = x;
s_bst = s;


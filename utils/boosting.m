function [x_bst, s_bst] = boosting(objfn, x_0, termfn, stepsize)
% Ivar's interpretation of a (perhaps not truly) boosting algorithm.
%
% ARGUMENTS: 
%    objfn     Objective function. Must accept vector x and return a scalar
%    x_0       Starting point for the optimization of vector x.
%    termfn    Termination function. Accepts arguments:
%                     n    Step number of this iteration.
%                     x    The present x being considered.
%                     s    Score as evaluated by objfn(x).
%    stepsize  Size of steps to take when optimizing x.
% 
% RETURNS:
%    x_bst     The best vector x found so far.
%    s_bst     The score of that vector, as evaluated by objfn.
%
% DETAILS:
% Given a starting vector x=x_0, this algorithm samples a step in each
% direction along every dimension in x and moves to that new point 
% if it improves the score according to objective function objfn(x) more
% than any alternative direction. 
%
% Will terminate whenever @termfn(n,x,s) returns a boolean true value, or
% the stepsize goes below 10^-9.
%
% Returns the best point found and the score of that point. 
%     
% Oct 12, 2012. Ivar Thorson.

% Ensure that the arguments are valid
[nr, nc] = size(x_0); 

if (nr ~= 1) error('x_0 must be a row vector'); end
if (stepsize <= 0) error('stepsize must be > 0'); end

% Starting search point
x = x_0;
s = objfn(x);
l = length(x_0);
n = 1;  % Step number

while stepsize > 10^-9 && ~termfn(n,x,s)
    x_pre = x;   % The state before taking any steps
    x_next = x;  % The best direction to step in so far
    s_next = s;
    
    % Try to take a step along every dimension
    fprintf('boosting.m: Stepping');
    
    for d = 1:l
        stepdir = zeros(1, l);
        stepdir(d) = stepsize;
                
        % Step in the direction
        x_fwd = x_pre + stepdir;
        x_bck = x_pre - stepdir;

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
    % Decrease the step size by 10x
     if all(x_pre == x)
         stepsize = stepsize * 0.1;
         fprintf('Stepsize decreased to %d\n', stepsize);
     else
        % Print the improvement after stepping
        fprintf('Scored: %d\n', s);
     end
    
    n = n + 1;
end

% Return the best value found so far
x_bst = x;
s_bst = s;


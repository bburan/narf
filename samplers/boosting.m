function [x_bst, s_bst] = boosting(x_0, objfn, termfn, stepsize)
% Ivar's interpretation of a (perhaps not truly) boosting algorithm.
%
% ARGUMENTS: 
%    x_0       Starting point for the optimization of vector x.
%    objfn     Objective function. Must accept vector x and return a scalar
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
% if it improves the score according to objective function objfn(x). 
%
% Will terminate whenever @termfn(n,x,s) returns a boolean true value.
%
% Start with as larger step size than you think you need, because the
% algorithm will decrease it as needed to improve fitting.
%
% Returns the best point found and the score of that point. 
%
% TODO:
%    Step size should probably be different for each dimension of x or
%    scaling problems could present difficulties. Or, maybe the variance
%    produced by each change should be the scaling factor?
%     
% Oct 12, 2012. Ivar Thorson.

% Ensure that the arguments are valid
[nr, nc] = size(x_0); 

if (nc ~= 1) error('x_0 must be a row vector'); end
if (stepsize <= 0) error('stepsize must be > 0'); end

% Starting search point
x = x_0;
s = objfn(x);
l = length(x_0);
n = 1;  % Step number

while ~termfn(n,x,s)
    x_pre = x;   % The state before taking any steps
    
    % Try to take a step along every dimension
    for d = 1:l
        stepdir = zeros(1, l);
        stepdir(d) = stepsize;
                
        x_fwd = x + stepdir;
        x_bck = x - stepdir;

        s_fwd = evalfn(x_fwd);
        s_bck = evalfn(x_bck);
        
        % Take a step forward if that is better
        if s_fwd < s
            s = s_fwd;
            x = x_fwd;
        end

        % Take a step backward if that is better
        if s_bck < s
            s = s_bck;
            x = x_bck;
        end
        
        fprintf('.');
    end
    
    % If the search point has not changed, step size is too big.
    % Decrease the step size by 10x
    if all(x_pre == x)
        stepsize = stepsize * 0.1;
        fprintf('log: stepsize = %d\n', stepsize);
    end
    
    % Print the improvement after stepping
    fprintf('log: new score: %d\n', s);

end

% Return the best value found so far
x_bst = x;
s_bst = s;


function [x_bst, s_bst, termcond] = boosting(objfn, x_0, termfn, stepsize, stepscale)
% [x_bst, s_bst, termcond] = boosting(objfn, x_0, termfn, stepsize=1, stepscale=10)
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

% Starting search point
x = x_0(:);
s = objfn(x);
l = length(x_0);
n = 0;  % Step number
s_delta = s;
while ~termfn(n, x, stepsize, s_delta)    
    x_next = x;  % x_next holds best stepping direction so far
    s_next = s;
    
    global XXX
    thismse=[XXX{end}.score_train_mse XXX{end}.score_test_mse];
    thiscorr=[XXX{end}.score_train_corr XXX{end}.score_test_corr];
    
    if n==0, mse1=[1 1]; mse0=thismse; corr1=thiscorr;  end
    mse1(n+2,:)=thismse./mse0;
    corr1(n+2,:)=thiscorr;
    fprintf('r_fit:    %8.4f      r_test:   %8.4f\n',thiscorr);
    fprintf('mse_fit:  %12.8f  mse_test:  %12.8f\n',thismse);
    fprintf('msen_fit: %12.8f  msen_test: %12.8f\n',...
            mse1(n+1,:)-mse1(n+2,:));
    if strcmp(XXX{1}.cellid,'oni009b-a1'),
        sfigure(1);clf;
        subplot(2,1,1);plot(mse1);legend('fit','test');
        title('normalized mse');
        subplot(2,1,2);plot(corr1);drawnow;
        title('corr coef');
    end
    
    
    % Try to take a step along every dimension
    fprintf('Boosting');
    
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
        % svd -- recalc error to reflect x_next in XXX
        s = objfn(x);
        
        % Print the improvement after stepping
        fprintf('coef# %3d, delta: %d, score:%d\n', dir, s_delta, s);
        n = n + 1;
        
     end
    
end

% Return the best value found so far
x_bst = x;
s_bst = s;
termcond = termfn(n, x, stepsize, s_delta);

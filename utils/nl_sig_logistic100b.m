function y = nl_sig_logistic100b(phi, z)
% The logistic sigmoid function, modified to allow asymmetetry.  
%  
% SVD ripped off of sig_logistic but curvature divided by 100 to make
% more tolerant to big step sizes
%
% JL changed the parameterization of this function to remove local minima
% induced by min/max and if conditions
    
    % Generic 5 parameter sigmoid parameters
    baserate = phi(1);
    peakrate = phi(2);
    if baserate > peakrate
        peakrate = baserate;
    end
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    lo = phi(4)/100;  % low side curvature
    if lo < 0
        lo=0;
    end
    hi = phi(5)/100;   % high side curvature
    if hi < 0
        hi = 0;
    end
    
    % Shift z left to right immediately
    z = z - lrshift;

    f = @(A, x) 1 ./ ( 1 + exp(-x*A));
    
    df = @(A, x) A*exp(-x*A) / (exp(-x*A)+1).^2;           

    opts = optimset('Display', 'off');
    
    if hi > lo
        q = (df(hi, 0.1) - df(lo, 0.1));
        if ~isfinite(q)
            inflection = Inf;
            ef = 0;
        else
            [inflection, ~, ef, ~] = fzero(@(x) df(hi,x) - df(lo,x), 0.1, opts);
        end
    else
        q = (df(hi, -0.1) - df(lo, -0.1));
        if ~isfinite(q)
            inflection = Inf;
            ef = 0;
        else
            [inflection, ~, ef, ~] = fzero(@(x) df(hi,x) - df(lo,x), -0.1, opts);
        end
    end
     
    if ef < 0
        % We reach this point typically when Phi is such that the sigmoid
        % has invalid parameters (such as when baserate=peakrate) and the
        % resulting sigmoid zeros everything to the same level. We can
        % safely ignore this event, since the fitters will realize that it
        % is a stupid step. 
        % keyboard;
    end
    
    offset = f(lo, inflection) - f(hi, inflection) ;

    y = z;
    ii = z<inflection;
    y(ii) = f(lo, z(ii));
        
    ii = z>=inflection;    
    y(ii) = f(hi, z(ii)) + offset;
   
    % Reset y to be between 0 and 1, since it got messed up with that offset
    y = y ./ (1 + offset);
    
    % Now scale and shift so baserate and peakrate are right
    y = (peakrate-baserate) * y + baserate;

end
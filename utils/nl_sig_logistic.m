function y = nl_sig_logistic(phi, z)
% The logistic sigmoid function, modified to allow asymmetetry

    % Generic 5 parameter sigmoid parameters
    baserate = min([phi(1), phi(2)]);
    peakrate = max([phi(1), phi(2)]);   
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    lo = phi(4);   % low side curvature
    hi = phi(5);   % high side curvature
    
    % Shift z left to right immediately
    z = z - lrshift;

    f = @(A, x) 1 ./ ( 1 + exp(-x*A));
    
    df = @(A, x) A*exp(-x*A) / (exp(-x*A)+1).^2;       
    
    if hi > lo
        inflection = fzero(@(x) df(hi,x) - df(lo,x), 0.1);
    else
        inflection = fzero(@(x) df(hi,x) - df(lo,x), -0.1);
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
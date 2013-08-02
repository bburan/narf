function y = nl_sig_cauchy(phi, z)
% The "Cauchy" sigmoid function, modified to allow asymmetry

    % Generic 5 parameter sigmoid parameters
    baserate = min([phi(1), phi(2)]);
    peakrate = max([phi(1), phi(2)]);   
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    lo = phi(4);   % low side curvature
    hi = phi(5);   % high side curvature term   
    
    % Shift z left to right immediately
    z = z - lrshift;    
    
    f = @(A, x) 2/pi * atan((pi/2) * A* x);
    
    df = @(A, x) A / ((pi^2 * A^2 * x.^2)*0.25 + 1);
    
    if hi > lo
        inflection = (2/pi) * -sqrt(1/(lo*hi));        
    else
        inflection = (2/pi) * -sqrt(1/(lo*hi));    
    end

    offset = f(lo, inflection) - f(hi, inflection) ;

    y = z;
    ii = z<inflection;
    y(ii) = f(lo, z(ii));
        
    ii = z>=inflection;    
    y(ii) = f(hi, z(ii)) + offset;
   
    % Reset y to be between -1 and 1, since it got messed up with that offset
    y = 2*((y + 1) ./ (2 + offset)) - 1;  
    
    % Now scale and shift so baserate, peakrate, and middle are fine
    y = (peakrate-baserate) * (0.5 * y + 0.5) + baserate;

end
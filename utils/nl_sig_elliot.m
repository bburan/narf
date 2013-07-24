function y = nl_sig_elliot(phi, z)
% The "Elliot" sigmoid function, modified to allow asymmetry
% From http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.46.7204
% "A better Activation Function for Artificial Neural Networks"

    % Generic 5 parameter sigmoid parameters
    baserate = min([phi(1), phi(2)]);
    peakrate = max([phi(1), phi(2)]);   
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    lo = phi(4);   % low side curvature
    hi = phi(5);   % high side curvature term   
    
    % Shift z left to right immediately
    z = z - lrshift;    
    
    if hi > lo
        inflection = -sqrt(1/(lo*hi));        
    else
        inflection = +sqrt(1/(lo*hi));       
    end    
    
    offset = lo*inflection / (1 + abs(lo*inflection)) - hi*inflection / (1 + abs(hi*inflection));
    
    y = z;
    ii = z<inflection;
    y(ii) = (lo*z(ii) ./ (1 + abs(lo*z(ii))));
        
    ii = z>=inflection;    
    y(ii) = (hi*z(ii) ./ (1 + abs(hi*z(ii)))) + offset;
   
    % Reset y to be between -1 and 1, since it got messed up with that offset
    y = 2*((y + 1) ./ (2 + offset)) - 1;  
    
    % Now scale and shift so baserate, peakrate, and middle are fine
    y = (peakrate-baserate) * (0.5 * y + 0.5) + baserate;

end
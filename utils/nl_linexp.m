function ret = nl_linexp(phi, z)
    baserate = phi(1);
    maxrate = phi(2);
    stimeffect = phi(3);
    if length(phi) > 3
        zthresh = phi(4);
    else
        zthresh = 0;
    end
    
    ret = baserate * ones(size(z));
    
    ii = z - zthresh > 0;    
    ret(ii) = maxrate - (maxrate-baserate)*exp(-abs(stimeffect) .* abs(z(ii) - zthresh));
    
end
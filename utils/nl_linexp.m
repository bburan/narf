function ret = nl_linexp(phi, z)
    baserate = phi(1);
    maxrate = phi(2);
    stimeffect = phi(3);    
    ret = maxrate - (maxrate-baserate)*exp(stimeffect .* abs(z)); 
    %ret = 1./(phi(1)^2 + phi(2) .* abs(z) + phi(3)); 
end
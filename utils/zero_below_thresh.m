function ret = zero_below_thresh(phi, z)   
    ret = (z - phi(1)) .* heaviside(z - phi(1));   
end
function ret = exponential(phi, z)
    ret = exp(phi(1) .* (z - phi(2))); 
end
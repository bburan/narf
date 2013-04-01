function ret = nl_zerothresh(phi, z)   
    ret = (z - phi(1)) .* heaviside(z - phi(1));   
end
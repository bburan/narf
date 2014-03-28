function ret = nl_exponential(phi, z)
    if length(phi) > 2
        offset = phi(3);
    else
        offset = 0;
    end
    ret = offset + exp(phi(1) .* (z - phi(2)));
end
function ret = nl_root(phi, z)
    exponent = phi(1);
    if length(phi) > 1
        baseline = phi(2);
    else
        baseline = 0;
    end
    ret = abs(z.^(1/exponent)) + baseline;
end

function ret = nl_root(phi, z)
    exponent = phi(1);
    if length(phi) > 1
        z0 = phi(2);
    else
        z0 = 0;
    end
    
    idxs = (z <= z0);
    z(idxs) = 0;
    ret = (z-z0).^(1/exponent);
    ret(idxs) = 0;
end

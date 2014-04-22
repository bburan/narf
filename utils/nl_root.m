function ret = nl_root(phi, z)
    exponent = phi(1);    
    if length(phi) > 1
        base = phi(2);
    else
        base = 0;
    end        
    if length(phi) > 2
        z0 = phi(3);
    else
        z0 = 0;
    end
    idxs = (z <= z0);
    z(idxs) = 0;
    ret = (z-z0).^(1/exponent) + base;
    ret(idxs) = base;
    
end

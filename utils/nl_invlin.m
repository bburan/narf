function ret = nl_invlin(phi, z)
    if length(phi) > 2
       thresh = phi(3); 
    else
       thresh = 0;
    end
    idx = z < thresh;
    ret = phi(1) + phi(2).*(1 - 1./((z-thresh) + 1)); 
    ret(idx) = phi(1);
end
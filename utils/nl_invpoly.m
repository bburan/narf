function ret = nl_invpoly(phi, z)
    ret = 1./(phi(1) .* z.^2 + phi(2) .* abs(z) + phi(3)); 
end
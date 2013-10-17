function ret = nl_zerothresh(phi, z)
    if length(phi)==1,
        ret = (z - phi(1)) .* heaviside(z - phi(1));   
    elseif length(phi)==2,
        ret = (z - phi(1)) .* heaviside(z - phi(1)).*phi(2);   
    elseif length(phi)==3,
        ret = (z - phi(1)) .* heaviside(z - phi(1)).*phi(2) +phi(3);   
    end
end
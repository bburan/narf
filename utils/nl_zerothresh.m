%function ret = nl_zerothresh(phi, z)
%
% phi [thresh gain offset] .  ph(2:3) optional
function ret = nl_zerothresh(phi, z)
    if length(phi)==1,
        ret = (z - phi(1)) .* heaviside(z - phi(1));   
    elseif length(phi)==2,
        ret = (z - phi(1)) .* heaviside(z - phi(1)).*phi(2)./100;   
    elseif length(phi)==3,
        ret = (z - phi(1)) .* heaviside(z - phi(1)).*phi(2)./100 +phi(3);   
    end
end
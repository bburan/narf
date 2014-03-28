%function ret = nl_zerothresh(phi, z)
%
% phi [thresh gain offset] .  ph(2:3) optional
% Ivar: Removed scaling terms because scaling by 1/100 results in too-large
% parameter values unless post-nonlinearity normalization is applied. 
function ret = nl_zerothresh(phi, z)
    if length(phi)==1,
        ret = (z - phi(1));
        ret(ret<0)=0;
    elseif length(phi)==2,
        ret = (z - phi(1));
        ret(ret<0)=0;
        ret=ret.*phi(2);
    elseif length(phi)==3,
        ret = (z - phi(1));
        ret(ret<0)=0;
        ret=ret.*phi(2)+phi(3);
    end
end
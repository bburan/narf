function ret = nl_log(phi, z)
% ret = nl_log(phi, z)
% phi(1) Curvature of the logarithm
% phi(2) Zero offset (base rate)
% phi(3) Input is ignored below this amount. 

    % The "offset" term actually determines the logarithm base. 
    offset = phi(1);
       
    % soften effects of more extreme offsets
    if offset>4,
        adjoffset=4+(offset-4)./10;
    elseif offset<-4
        adjoffset=-4+(offset+4)./10;
    else
        adjoffset=offset;
    end
    
    % Offset from zero
    if length(phi) > 1
        zeroer = phi(2);
    else
        zeroer = 0;
    end   
    
    % Zero below threshold
    if length(phi) > 2
        zbt = phi(3);
        z(z<zbt) = zbt;
        z = z-zbt;
    end
    
    ret = log((z + 10^adjoffset)/10^adjoffset) + zeroer;
end

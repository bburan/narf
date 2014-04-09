function ret = nl_dlog(phi, z)
% ret = nl_log(phi, z)
% phi(1) Curvature of the logarithm
% phi(2) Zero offset (base rate)
% phi(3) Input is ignored below this amount. 

    offset = phi(1);
    
    % soften effects of more extreme offsets
    if offset>4,
        adjoffset=4+(offset-4)./10;
    elseif offset<-4
        adjoffset=-4+(offset+4)./10;
    else
        adjoffset=offset;
    end
    
    d = 10^adjoffset;    
    
    % Offset from zero
    if length(phi) > 1
        zeroer = phi(2);
    else
        zeroer = -log(d);
    end       
    
    % Zero below threshold
    if length(phi) > 2
        zbt = phi(3);
    else        
        zbt = 0;
    end      
    
    % If any values are below zero, make them 0    
    z(z<zbt) = zbt;
    z = z-zbt;   
    
    ret = log(log(z+d)-log(d)+d) + zeroer;
end

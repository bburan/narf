function ret = nl_dlog(phi, z)
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
    
    ret = log(log(z+d)-log(d)+d);
end

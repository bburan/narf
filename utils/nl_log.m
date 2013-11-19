function ret = nl_log(phi, z)
    offset = phi(1);
    if length(phi) > 1
        zeroer = phi(2);
    else
        zeroer = -log(10^offset);
    end
    
    % soften effects of more extreme offsets
    if offset>4,
        adjoffset=4+(offset-4)./10;
    elseif offset<-4
        adjoffset=-4+(offset+4)./10;
    else
        adjoffset=offset;
    end
    
    ret = log(z + 10^adjoffset) + zeroer;
end

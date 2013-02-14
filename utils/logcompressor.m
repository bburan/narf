function ret = logcompressor(phi, z)
    offset = phi(1);
    if length(phi) > 1
        zeroer = phi(2);
    else
        zeroer = 0;
    end
    ret = log(z + 10^offset) + zeroer;
end

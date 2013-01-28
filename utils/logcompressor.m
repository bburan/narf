function ret = logcompressor(phi, z)
    offset = phi(1);
    ret = log(z + 10^offset);
end

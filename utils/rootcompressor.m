function ret = rootcompressor(phi, z)
    exponent = phi(1);
    ret = abs(z.^(1/exponent));
end

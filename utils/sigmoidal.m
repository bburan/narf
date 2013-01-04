function ret = sigmoidal(phi, z)
    mu = phi(1);
    sigma = phi(2);
    amp = phi(3);
    offset = phi(4);
    
    ret = amp * (1 + erf((z - mu) / sqrt(2 * sigma^2))) + offset;
end

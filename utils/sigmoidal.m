function ret = sigmoidal(phi, z)
    mu = phi(1);
    sigma = phi(2);
    amp = phi(3);
    
    ret = 0.5 * (1 + erf((z - mu) / sqrt(2 * sigma^2)));
end
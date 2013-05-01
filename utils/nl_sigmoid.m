function ret = nl_sigmoid(phi, z)
    mu = phi(1);
    sigma = phi(2);
    amp = phi(3);
    offset = phi(4);
    % The abs() makes it so that all sigmoids are increasing from L to R
    %ret = abs(amp) * (1 + erf((z - mu) / sqrt(2 * sigma^2))) + offset;
    ret = amp * (1 + erf((z - mu) / sqrt(2 * sigma^2))) + offset;
end

function y = nl_richards(phi, z)
% The richards growth equation sigmoid

    base = min([phi(1), phi(2)]); 
    peak = max([phi(1), phi(2)]);
    gamma = phi(3);   % Where the centerpoint of the sigmoid is
    kappa = phi(4);   % Curvature
    bias  = phi(5);   % The left/right bias. I actually take the exponent so it can be positive or negative, and 0 is no bias.
       
    y = base + (peak-base) ./ (1 + exp(bias).*exp(-kappa.*(z-gamma))).^(1/(exp(bias)));
end



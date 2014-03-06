function y = nl_dexp(phi, z)
% The double exponential sigmoid, naturally asymmetric (Aka gompertz?)

    base = min([phi(1), phi(2)]); 
    peak = max([phi(1), phi(2)]);
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    kappa = phi(4);   % Curvature
    
    y = peak * exp(-exp(-kappa.*(z-lrshift))); 

end
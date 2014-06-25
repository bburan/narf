function y = nl_dexp2(phi, z)
% The double exponential sigmoid with a different parameterization
% (the min/max is not a very good idea for fitting, as it effectively
% creates 2 global minima -- something that we want to avoid at all cost
% becaues it makes the exploration of the parameter space more difficult)

    base = phi(1);
    peak = phi(2);
    if peak < 0
        % we don't allow "reverse" sigmoid
        y = base;
        return
    end
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    kappa = phi(4);   % Curvature
    
    y = base + peak * exp(-exp(-kappa.*(z-lrshift))); 
    
end
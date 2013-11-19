function ret = nl_softzero(phi, z)
    alpha = phi(1);  % Scaling parameter
    beta  = phi(2);  % Curvature
    theta = phi(3);  % Offset   
    
    if (length(phi) > 3)
        yoffset = phi(4);  % Offset   
    else
        yoffset = 0;
    end
    ret = alpha * log(1 + exp(beta*(z-theta))) + yoffset; 
end

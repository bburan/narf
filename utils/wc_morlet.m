function weights = wc_morlet(phi, N_inputs)
    [N_chans, N_parms] = size(phi);
    
    if N_parms ~= 3
        error('WC_MORLET needs exactly three parameters per channel');
    end
    
    mu_khz  = phi(:, 1);  % Center of gaussian in kHz
    a       = phi(:, 2);  % Scaling factor
    k0      = phi(:, 3);  % ??  
    
    weights = zeros(N_inputs, N_chans);
    
    for c = 1:N_chans
        if (mu_khz(c) < 0.2 || mu_khz(c) > 20)
            % If you go below 200Hz or above 20000Hz, make the prediction terrible
            weights(:, c ) = ones(N_inputs, 1);
        else
            % Otherwise, position the gaussian in the right place
            spacing = (log(20000) - log(200)) / N_inputs;
            mu = (log(mu_khz(c)*1000) - log(200)) / spacing;            
            x = 1:N_inputs;
            weights(:, c) = real(exp(-(x-mu).^2 / a(c)^2) .* exp(-1i * k0(c) * (x-mu)));
        end
    end
end

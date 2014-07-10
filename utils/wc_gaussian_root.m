function weights = wc_gaussian_root(phi, N_inputs)
    [N_chans, N_parms] = size(phi);
    
    if N_parms ~= 3 
        error('WC_GAUSSIAN needs exactly 3 parameters per channel');
    end
    
    mu_khz  = phi(:, 1);  % Center of gaussian in kHz
    sig_khz = phi(:, 2);  % Gaussian STDDEV (as if kHz were linear)
    nroot   = phi(:, 3);  % Nth root 
    
    weights = zeros(N_inputs, N_chans);
    
    for c = 1:N_chans
        if (mu_khz(c) < 0.2 || mu_khz(c) > 20)
            % If you go below 200Hz or above 20000Hz, make the prediction terrible
            weights(:, c ) = ones(N_inputs, 1);
        else
            % Otherwise, poisition the gaussian in the right place
            spacing = (log(20000) - log(200)) / N_inputs;
            mu = (log(mu_khz(c)*1000) - log(200)) / spacing;
            sigma = (sig_khz(c)/10) * mu;
            weights(:, c) = gauss1([mu, sigma], 1:N_inputs).^(1/nroot(c));
        end
    end
end

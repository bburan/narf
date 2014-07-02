function y = nl_sig_5PL(phi, z)
% The reparameterized 5 PL (5 parameter-logistic) sigmoid function, used
% frequently in dose-response curves. This is an extension of the 4 PL 
% which allows assymetric curves.
% (see http://onlinelibrary.wiley.com/doi/10.1002/cem.1218/pdf )
%   Parameters:
%   => phi(1) is the baserate and should be positive
%   => phi(2) is the amplitude (maxrate-baserate) and should be positive
%   => phi(3) is the input that triggers the mean firing rate and should be
%   positive
%   => phi(4) is the curvature parameter and should be positive (values in [0.01 100] make sense)
%   => phi(5) is the assymetry parameter and should be positive (values in [0.01 100] make sense)

    % Generic 5 parameter sigmoid parameters
    baserate = phi(1);
    peakrate = phi(1)+phi(2);
    C = phi(3);
    curvature = phi(4);
    assymetry = phi(5);
%     assymetry = 1;
    
    z(z<0) = 0;
    
    if baserate < 0, baserate = 0; end;
    if peakrate < baserate, peakrate = baserate; end;
    if C < 0, y=0*z; return; end;
    if curvature < 0, y=0*z; return; end;
    if assymetry < 0, y=0*z; return; end;
    
    if baserate > 100, baserate = 100; end;
    if peakrate > 100, peakrate = 100; end;
    if C >100, C=100; end;
    if curvature > 100, curvature=100; end;
    if assymetry > 100, assymetry=100; end;

    y = peakrate - (peakrate-baserate) ./ (1 + (2^(1/assymetry) -1 ) * ( (z/C).^curvature ) ).^assymetry;

end
function y = nl_sig_gumber(phi, z)
% The left gumbel sigmoid function, which is asymmetric 

    % Generic 5 parameter sigmoid parameters
    baserate = min([phi(1), phi(2)]);
    peakrate = max([phi(1), phi(2)]);   
    lrshift = phi(3); % Where the centerpoint of the sigmoid is
    mu = phi(4);   % low side curvature
    beta = phi(5);   % high side curvature      

    % Shift z left to right immediately
    % The median value of a gumbel distribution is 
    %     mu - beta*ln(ln(2)),
    % so shifting is a little more complicated than normal
    z = z - lrshift + (mu-beta*log(log(2)));
   
    f = @(mu, beta, x) 1 - exp(-exp((x-mu)./beta));
    
    y = f(mu, beta, z);
   
    % Now scale and shift so baserate and peakrate are right
    y = (peakrate-baserate) * y + baserate;

end
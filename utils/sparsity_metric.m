function sparsity = sparsity_metric(coefs)
    % A metric of the inverse sparsity (perhaps this should be called the 
    % relative dimensionality or something? 
    % Varies between 1 (all coefs) and N (the number of coefficients)
    % Thus, it is a comparable metric regardless of how many coefficients
    % there are: if the number is 3.1, then it has a dimensionality or 
    % complexity of about 3 strong coefficients (and everything else is
    % assumed to be nearly zero).

    ncoefs = numel(coefs);
    
    if ncoefs < 2
        sparsity = 1;
        return
    end
  
    c_bar = sqrt(sum(coefs(:).^2)); % L2 norm
    sparsity = sum(abs(coefs(:)) / c_bar)^2;

end
% 
% % Check the bounds of the function with this
% clear
%  for ii = 1:100,
%      X(ii) = ii;
%      nosignal = ones(1, ii);
%      Y1(ii) = sparsity_metric(nosignal);
%      perfect = zeros(1, ii);
%      perfect(1) = 1;
%      Y2(ii) = sparsity_metric(perfect);
%  end
% 
% plot(X, Y1, X, Y2);
% legend('No signal', 'Perfect');
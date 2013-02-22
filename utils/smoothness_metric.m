function smoothness = smoothness_metric(coefs)
    smoothness = filter([1,-1], 1, coefs, [], 2);
    smoothness = sqrt(sum(smoothness(:).^2));
end
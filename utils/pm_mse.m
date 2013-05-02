function score = pm_mse()
% pm_mse()
%
% Performance metric: MSE with optional sparseness penalty

global XXX META;

if ~isfield(META, 'sparsity_weight')
    score = XXX{end}.score_train_mse; 
else
    [fir_mods, fir_idxs] = find_modules('fir_filter');
    sparsities = zeros(size(fir_idxs));   
    for ii=1:length(fir_idxs)
        sparsities(ii) = sparsity_metric(fir_mods{ii}.coefs);
    end
    score = XXX{end}.score_train_mse + META.sparsity_weight * sum(sparsities);
end
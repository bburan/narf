function score = pm_mse()
% pm_mse()
%
% Performance metric: MSE with optional sparseness penalty

global XXX META STACK;

if ~isfield(META, 'sparsity_weight')
    score = XXX{end}.score_train_mse; 
else
    [fir_mods, fir_idxs] = find_modules(STACK, 'fir_filter');
    sparsities = [];
    for ii=1:length(fir_idxs)
        for pp=1:length(fir_mods{ii})
            sparsities(end+1) = sparsity_metric(fir_mods{ii}{pp}.coefs);
        end
    end
    score = XXX{end}.score_train_mse + (META.sparsity_weight * sum(sparsities));
end
function [mse, penalty] = pm_mse()
% pm_mse()
%
% Performance metric: MSE with optional sparseness penalty
%
% RETURNS:
%    mse      Mean Squared Error
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

if ~isfield(META, 'sparsity_weight')
    mse = XXX{end}.score_train_mse; 
else    
    [fir_mods, fir_idxs] = find_modules(STACK, 'fir_filter');
    sparsities = [];
    for ii=1:length(fir_idxs)
        for pp=1:length(fir_mods{ii})
            sparsities(end+1) = sparsity_metric(fir_mods{ii}{pp}.coefs);
        end
    end
    mse = XXX{end}.score_train_mse;
    penalty = (META.sparsity_weight * sum(sparsities));
end
function [mse, penalty, val_mse] = pm_mse()
% pm_mse()
%
% Performance metric: MSE with optional sparseness penalty
%
% RETURNS:
%    mse      Mean Squared Error
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

mse = XXX{end}.score_train_mse; 
val_mse = XXX{end}.score_test_mse; 

if ~isfield(META, 'sparsity_weight')
    penalty = 0;
else    
    [fir_mods, fir_idxs] = find_modules(STACK, 'fir_filter');
    sparsities = [];
    for ii=1:length(fir_idxs)
        for pp=1:length(fir_mods{ii})
            sparsities(end+1) = sparsity_metric(fir_mods{ii}{pp}.coefs);
        end
    end
    penalty = (META.sparsity_weight * sum(sparsities));
end
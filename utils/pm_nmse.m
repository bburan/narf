function [nmse, penalty, val_nmse] = pm_nmse()
% pm_nmse()
%
% Performance metric: Normalized MSE with optional sparseness penalty
%
% RETURNS:
%    nmse      Mean Squared Error
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

nmse = XXX{end}.score_train_nmse; 
val_nmse = XXX{end}.score_test_nmse; 

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
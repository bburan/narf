function [norm, penalty, val_norm] = pm_norm()
% pm_norm()
%
% Performance metric: MSE with optional sparseness penalty
%
% RETURNS:
%    norm     P-Norm
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

norm = XXX{end}.score_train_norm; 
val_norm = XXX{end}.score_test_norm;

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
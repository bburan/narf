function [norm, penalty] = pm_norm()
% pm_norm()
%
% Performance metric: MSE with optional sparseness penalty
%
% RETURNS:
%    norm     P-Norm
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

if ~isfield(META, 'sparsity_weight')
    norm = XXX{end}.score_train_norm; 
    penalty = 0;
else    
    [fir_mods, fir_idxs] = find_modules(STACK, 'fir_filter');
    sparsities = [];
    for ii=1:length(fir_idxs)
        for pp=1:length(fir_mods{ii})
            sparsities(end+1) = sparsity_metric(fir_mods{ii}{pp}.coefs);
        end
    end
    norm = XXX{end}.score_train_norm;
    penalty = (META.sparsity_weight * sum(sparsities));
end
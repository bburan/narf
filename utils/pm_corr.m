function [corr, penalty, val_corr] = pm_corr()
% pm_corr()
%
% Performance metric: Correlation with spareseness penalty
%
% RETURNS:
%    corr     One minus the Correlation
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

corr = 1 - XXX{end}.score_train_corr; 
val_corr = 1 - XXX{end}.score_val_corr;

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
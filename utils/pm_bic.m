function [bicinv, penalty, val_bicinv] = pm_bic()
% pm_bic()
%
% Performance metric: 
%
% RETURNS:
%    norm     P-Norm
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

bic = XXX{end}.score_train_bic; 
bicinv = 1/bic;
val_bic = XXX{end}.score_test_bic; 
val_bicinv = 1/val_bic; 

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
function [bicinv, penalty] = pm_bic()
% pm_bic()
%
% Performance metric: 
%
% RETURNS:
%    norm     P-Norm
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

if ~isfield(META, 'sparsity_weight')
    bic = XXX{end}.score_train_bic; 
    penalty = 0;
    bicinv = 1/bic;
else    
    [fir_mods, fir_idxs] = find_modules(STACK, 'fir_filter');
    sparsities = [];
    for ii=1:length(fir_idxs)
        for pp=1:length(fir_mods{ii})
            sparsities(end+1) = sparsity_metric(fir_mods{ii}{pp}.coefs);
        end
    end
    bic = XXX{end}.score_train_bic;
    penalty = (META.sparsity_weight * sum(sparsities));
    bicinv = 1/bic;
end
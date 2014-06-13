function [nlogl, penalty, val_nlogl] = pm_nlogl()
% pm_nlogl()
%
% Performance metric: 
%
% RETURNS:
%    nlogl
%    penalty  Sparsity of the solution multiplied by META.sparsity_weight

global XXX META STACK;

nlogl = XXX{end}.score_train_nlogl;
val_nlogl = XXX{end}.score_test_nlogl; 

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
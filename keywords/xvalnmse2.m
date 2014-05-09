function xvalnmse2(crossvalidation_fold)

global MODULES META STACK;


if ~exist('crossvalidation_fold', 'var')
    crossvalidation_fold = 2;
end

% mods = find_modules(STACK, 'mean_squared_error', true);
% if isempty(mods)
%     append_module(MODULES.mean_squared_error.mdl(struct('crossvalidation_fold', crossvalidation_fold)));   
% end

mods = find_modules(STACK, 'passthru', true);
if isempty(mods)
    append_module(MODULES.passthru.mdl(struct('crossvalidation_fold', crossvalidation_fold)));   
else
    if isfield(mods,'crossvalidation_fold'),
        mods.crossvalidation_fold = crossvalidation_fold;
    end
end

META.perf_metric = @pm_nmse;
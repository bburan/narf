function mixfit3firperfile()

global STACK

mixfit3(); %% Relative boost algorithm across all files to initialize

[~, mod_idxs] = find_modules(STACK, 'fir_filter', false);
for ii=1:(mod_idxs{1}(1)-1),
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields={};
    end
end

% Then boost on each file individually
fit_split_simply(@mixfit3nomse, @split_by_respfile, @unify_respfiles); 



function mixfit3nlpercode()

global STACK

% We need a performance metric, so add one
mse();

% Relative boost algorithm across all files to initialize
mixfit3nomse();

[~, mod_idxs] = find_modules(STACK, 'fir_filter', false);
for ii=1:mod_idxs{end}(end),
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields={};
    end
end

% Then boost on each file individually
fit_split_simply(@mixfit3nomse, @split_by_filecode, @unify_respfiles); 



function modulekeys = module_block_combos(mm)
% ARGUMENTS:
%   MM   A module blocks structure like the first argument of fit_models()
%
% RETURNS:
%   MODULEKEYS        A cell array of modulekeys

opts = cellfun(@fieldnames, mm, 'UniformOutput', false);
N_opts = cellfun(@(m) length(fieldnames(m)), mm);
N_models = prod(N_opts);
modulekeys     = cell(N_models,1);

for ii = 1:N_models,
    % Behold matlab's ugliness for destructuring binds!
    [i1 i2 i3 i4 i5 i6 i7 i8 i9] = ind2sub(N_opts, ii);
    indexes = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    
    % Build the model name
    opt_names = arrayfun(@(gi, ind) opts{gi}{ind}, ...
        1:length(opts), indexes(1:length(opts)), ...
        'UniformOutput', false);
    
    modulekeys{ii} = opt_names;

end
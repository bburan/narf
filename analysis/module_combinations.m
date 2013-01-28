function [models, modelnames] = module_combinations(mm)
% ARGUMENTS:
%   MM   A module blocks structure like the first argument of fit_models()
%
% RETURNS:
%   MODELS        A cell array of models
%   MODELNAMES    A cell array of model names

opts = cellfun(@fieldnames, mm, 'UniformOutput', false);
N_opts = cellfun(@(m) length(fieldnames(m)), mm);
N_models = prod(N_opts);
models     = cell(N_models,1);
modelnames = cell(N_models,1);
for ii = 1:N_models,
    % Behold matlab's ugliness for destructuring binds!
    [i1 i2 i3 i4 i5 i6 i7 i8 i9] = ind2sub(N_opts, ii);
    indexes = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    
    % Build the model name
    opt_names = arrayfun(@(gi, ind) opts{gi}{ind}, ...
        1:length(opts), indexes(1:length(opts)), ...
        'UniformOutput', false);
    tmp = cellfun(@(n) sprintf('%s_', n), opt_names, 'UniformOutput', false);
    
    modelname = strcat(tmp{:});
    modelnames{ii} = modelname(1:end-1); % Remove trailing underscore
    
    % Build the model 
    blocks = cellfun(@getfield, mm, opt_names, 'UniformOutput', false);
    models{ii} = {};
    for jj = 1:length(blocks),
        for kk = 1:length(blocks{jj})
            models{ii}{end+1} = blocks{jj}{kk};
        end
    end  
    
end
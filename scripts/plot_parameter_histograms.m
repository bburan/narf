function plot_parameter_histograms(batch, cellids, modelnames)

stack_extractor = @pack_fittables;

[params, ~, ~] = load_model_batch(batch, cellids, modelnames, ...       
                                  stack_extractor);

% Check that the number of free parameters is the same in all models
n_params = [];
for ii = 1:numel(params)
    prms = params{ii};
    if isempty(n_params) 
        n_params = length(prms);
    else
        if isempty(prms)
            continue;
        end
        if n_params ~= length(prms)
            error(['plot_parameter_histograms.m requires all models ' ...
                   'to have the same number of free parameters.']);
        end
    end
end

figure('Name', 'Parameter Histograms', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
n = ceil(sqrt(n_params));

for ii = 1:n_params
    subplot(n, n, ii);
    data = cellfun(@(c) c(ii), params);
    hist(data, 50);
    xlabel(sprintf('Parameter %d', ii));
end


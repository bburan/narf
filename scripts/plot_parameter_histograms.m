function plot_parameter_histograms(batch, cellids, modelnames)

% JL June 2014: added this function to get the names in addition to the
% values (mostly copy+paste from pack_fittable.m)
    function w_names = pack_fittables_and_names(~)
        global STACK;
        w = [];
        names = {};
        w_names = cell(2,1);
        for iii = 1:length(STACK)
            mm = STACK{iii};
            nsplits = length(mm);
            for kk = 1:nsplits
                m = mm{kk};
                if isfield(m, 'fit_fields')
                    for jj = 1:length(m.fit_fields),
                        p = m.fit_fields{jj};
                        for ll = 1:numel(m.(p))
                            names = {names{:} [p sprintf('[%d]',ll)]};
                        end
                        w = cat(1, w, reshape(m.(p), numel(m.(p)), 1));
                    end
                end
            end
        end
        w_names{1} = w;
        w_names{2} = names;
    end

% stack_extractor = @pack_fittables;
stack_extractor = @pack_fittables_and_names;


[params, ~, ~] = load_model_batch(batch, cellids, modelnames, ...
    stack_extractor);

% Check that the number of free parameters is the same in all models
n_params = [];
for ii = 1:numel(params{1})
    prms = params{1}{ii};
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

for ii = numel(params):-1:1
    if isempty(params{ii})
        params(ii) = [];
        fprintf('excluding empty model %d\n',ii)
    end
end


figure('Name', ['Parameter Histograms of ' cell2mat(modelnames(1))], 'NumberTitle', 'off', 'Position', [20 50 900 900]);
n = ceil(sqrt(n_params));

warning('off','all');
try
    dat = cellfun(@(c) mat2cell(c{1}), params);
    csvwrite('/tmp/parameter_histogram.csv', dat);
catch whatever
    fprintf('could not coexerce data to export to csv (cause: %s)\n', getReport(whatever));
end
warning('on','all');

for ii = 1:n_params
    subplot(n, n, ii);
    data = cellfun(@(c) c{1}(ii), params);
    unique_names = char(unique(cellfun(@(c) c{2}(ii), params))');
    hist(data, 50);
    xlabel(sprintf('%s [%3.2f to %3.2f]', unique_names, min(data), max(data)));
end

end

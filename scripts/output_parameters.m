 function output_parameters(batch, cellids, modelnames)

% JL June 2014: added this function to get the names in addition to the
% values (mostly copy+paste from pack_fittable.m)
    function w_names = pack_useful_info(~)
        global STACK XXX;
        w = [];
        names = {};
        w_names = cell(2,1);
        K = 0;
        

        for iii = 1:length(STACK)
            mm = STACK{iii};
            nsplits = length(mm);
            for kk = 1:nsplits
                m = mm{kk};
                if isfield(m, 'fit_fields')
                    for jj = 1:length(m.fit_fields),
                        p = m.fit_fields{jj};
                        K = K + numel(m.(p));
                        for ll = 1:numel(m.(p))
                            names = {names{:} [p sprintf('_%d',ll)]};
                        end
                        w = cat(1, w, reshape(m.(p), numel(m.(p)), 1));
                    end
                end
            end
        end
        
        calc_xxx(1);
        
        n = 0;
        for tt = 1:length(XXX{end}.training_set)
            n = n + sum(sum(~isnan(XXX{end}.dat.(XXX{end}.training_set{tt}).respavg)));
        end
        
        w = [w; XXX{end}.score_train_nmse; XXX{end}.score_test_nmse; ...
            XXX{end}.score_train_corr; XXX{end}.score_test_corr; ...
            n * log(XXX{end}.score_train_mse) + 2*K; n * log(XXX{end}.score_train_mse) + 2*K + (2*K*(K+1))/(n-K-1); ...
            n * log(XXX{end}.score_test_mse) + 2*K; n * log(XXX{end}.score_test_mse) + 2*K + (2*K*(K+1))/(n-K-1) ];
        names = {names{:} 'train_score' 'test_score' 'train_corr' 'test_corr' 'train_AIC' 'train_AICC' 'test_AIC' 'test_AICC'};
        
        w_names{1} = w;
        w_names{2} = names;
    end

stack_extractor = @pack_useful_info;

for n_i = 1:numel(modelnames)
    
%     filename = '/tmp/parameter_histogram.csv';

    filename = ['/auto/users/lienard/data/params_', num2str(batch), '_', cell2mat(modelnames(n_i)), '.csv'];

%     [FileName,PathName,~] = ...
%         uiputfile(['/auto/users/lienard/data/params_', num2str(batch), '_', cell2mat(modelnames(n_i)), '.csv'],'Save to...');
%     if ~FileName,
%         disp('Canceled save.');
%         continue;
%     else
%         fprintf('Saving csv to %s%s\n',PathName,FileName);
%         filename = [PathName FileName];
%     end

    fid = fopen(filename, 'w');


    [params, ~, ~] = load_model_batch(batch, cellids, {modelnames{n_i}}, ...
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
    
    warning('off','all');
    dat = cellfun(@(c) mat2cell(c{1}), params);
    warning('on','all');
    
    unique_names = char(unique(cellfun(@(c) c{2}(1), params))');
    for ii = 2:n_params
        unique_names = [unique_names ',' char(unique(cellfun(@(c) c{2}(ii), params))')];
    end
    
    fprintf(fid, '%s\n', unique_names);
    fprintf('%s\n', unique_names);
    fclose(fid);
    dlmwrite(filename, cell2mat(dat)', '-append');
    
end

end

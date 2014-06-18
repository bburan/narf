function compute_akaike(batch, cellids, modelnames)
% Computes several interesting graphs for a combination of models

global XXX STACK;

if nargout==0,
    [FileName,PathName,~] = ...
        uiputfile(['akaike_', num2str(batch), '.csv'],'Save to...');
    if ~FileName,
        disp('Canceled save.');
        return;
    else
        fprintf('Saving csv to %s%s\n',PathName,FileName);
        fid = fopen([PathName FileName],'w');
        fprintf(fid,'CELL_ID MODEL TRAINING_SCORE TRAINING_CORR TRAINING_AIC TRAINING_AICC TESTING_SCORE TESTING_CORR TESTING_AIC TESTING_AICC\n');
    end
end

[~, metas] = load_model_batch(batch, cellids, modelnames);

n_models = size(metas,1);
n_cells = size(metas,2);

training_score = zeros(n_models,n_cells);
training_corr = zeros(n_models,n_cells);
training_AIC = zeros(n_models,n_cells);
training_AICc = zeros(n_models,n_cells);
testing_score = zeros(n_models,n_cells);
testing_corr = zeros(n_models,n_cells);
testing_AIC = zeros(n_models,n_cells);
testing_AICc = zeros(n_models,n_cells);


for c = 1:n_cells
    for m = 1:n_models
        
        try
            load_model(metas{m,c}.modelpath);
            
%             STACK{1}{1}.stimulus_channel_count = 36;
%             for ii=1:length(STACK)
%                 if isfield(STACK{ii}{1},'name')
%                     if strcmp(STACK{ii}{1}.name, 'lindeberg_filter')
%                         STACK{ii}{1}.lincoefs(1) = 2 * ...
%                             STACK{ii}{1}.lincoefs(1);
%                         STACK{ii}{1}.lincoefs(3) = 2 * ...
%                             STACK{ii}{1}.lincoefs(3);
%                         STACK{4}{1}.num_dims = 36;
%                     end
%                 end
%             end
            
            calc_xxx(1);
            
            crossval = 0;
            K = 0;
            n = 0;
            for tt = 1:length(XXX{end}.training_set)
                n = n + sum(sum(~isnan(XXX{end}.dat.(XXX{end}.training_set{tt}).respavg)));
                %         s = size(XXX{end}.dat.(XXX{end}.training_set{tt}).respavg);
                %         n = n + s(1) * s(2);
            end
            
            for ii=1:length(STACK)
                if isfield(STACK{ii}{1},'fit_fields')
                    fitfields=STACK{ii}{1}.fit_fields;
                    for jj=1:length(fitfields)
                        for kk=1:length(STACK{ii})
                            K = K + length(STACK{ii}{kk}.(fitfields{jj}));
                        end %kk
                    end %jj
                end %if
                if isfield(STACK{ii}{1},'name')
                    if strcmp(STACK{ii}{1}.name, 'lindeberg_filter')
                        % we add one meta-parameter for the lindeberg filter
                        K = K + 1;
                    end
                    if strcmp(STACK{ii}{1}.name, 'passthru')
                        % we handle differently the cross-validation case
                        if isfield(STACK{ii}{1},'crossvalidation_fold')
                            if STACK{ii}{1}.crossvalidation_fold ~= 0,
                                crossval = 1;
                            end
                        end
                    end
                end
            end
            
            if crossval,
                % for now, we must ignore the correlation as it is computed on
                % the whole training set and not on the right fold.
                training_corr(m,c) = 0;
                testing_corr(m,c) = 0;
                n = n / 2;
            else
                training_corr(m,c) = XXX{end}.score_train_corr;
                testing_corr(m,c) = XXX{end}.score_test_corr;
            end
            
            training_score(m,c) = XXX{end}.score_train_nmse;
            testing_score(m,c) = XXX{end}.score_test_nmse;
            
            training_AIC(m,c) = n * log(XXX{end}.score_train_mse) + 2*K;
            training_AICc(m,c) = n * log(XXX{end}.score_train_mse) + 2*K + (2*K*(K+1))/(n-K-1);
            testing_AIC(m,c) = n * log(XXX{end}.score_test_mse) + 2*K;
            testing_AICc(m,c) = n * log(XXX{end}.score_test_mse) + 2*K + (2*K*(K+1))/(n-K-1);
            
            fprintf(fid,'%s %s %f %f %f %f %f %f %f %f\n', cellids{c} , modelnames{m}, training_score(m,c), training_corr(m,c), training_AIC(m,c), training_AICc(m,c), testing_score(m,c), testing_corr(m,c), testing_AIC(m,c), testing_AICc(m,c));
            fprintf(1,'%s %s %f %f %f %f %f %f %f %f\n', cellids{c} , modelnames{m}, training_score(m,c), training_corr(m,c), training_AIC(m,c), training_AICc(m,c), testing_score(m,c), testing_corr(m,c), testing_AIC(m,c), testing_AICc(m,c));
            
        catch err
            fprintf(1,'catched an error\n');
        end
    end
end

% for c = 1:n_cells
%     for m = 1:n_models
%         fprintf(fid,'%s %s %f %f %f %f\n', cellids{c} , modelnames{m}, training_score(m,c), training_corr(m,c), AIC(m,c), AICc(m,c), testing_score(m,c), testing_corr(m,c));
%     end
% end

end


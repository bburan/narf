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
        fprintf(fid,'CELL_ID MODEL TRAINING_SCORE TRAINING_CORR AIC AICC\n');
    end
end

[~, metas] = load_model_batch(batch, cellids, modelnames);

n_models = size(metas,1);
n_cells = size(metas,2);

training_score = zeros(n_models,n_cells);
training_corr = zeros(n_models,n_cells);
AIC = zeros(n_models,n_cells);
AICc = zeros(n_models,n_cells);

for c = 1:n_cells
    for m = 1:n_models
        load_model(metas{m,c}.modelpath);
        calc_xxx(1);
        train_name = XXX{end}.training_set{1};
        s = size(XXX{end}.dat.(train_name).stim);
%         n = s(1) * s(2);
        n = s(2);
        K = 0;
        
        for ii=1:length(STACK)
            if isfield(STACK{ii}{1},'fit_fields')
                fitfields=STACK{ii}{1}.fit_fields;
                for jj=1:length(fitfields)
                    for kk=1:length(STACK{ii})
                        K = K + length(STACK{ii}{kk}.(fitfields{jj}));
                    end
                end
            end
        end
        
        training_score(m,c) = XXX{end}.score_train_nmse;
        training_corr(m,c) = XXX{end}.score_train_corr;
        AIC(m,c) = n * log(XXX{end}.score_train_nmse) + 2*K;
        AICc(m,c) = n * log(XXX{end}.score_train_nmse) + 2*K + (2*K*(K+1))/(n-K-1);
        
    end
end

for c = 1:n_cells
    for m = 1:n_models
%         fprintf(1,'AIC %s = %f\n', modelnames{m,c}, AIC(m,c));
        fprintf(fid,'%s %s %f %f %f %f\n', cellids{c} , modelnames{m}, training_score(m,c), training_corr(m,c), AIC(m,c), AICc(m,c));
    end
end

end


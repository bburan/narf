function compute_akaike(batch, cellids, modelnames)
% Computes several interesting graphs for a combination of models

global XXX STACK;

[~, metas] = load_model_batch(batch, cellids, modelnames);

n_models = size(metas,1);
n_cells = size(metas,2);

AIC = zeros(n_models,n_cells);
AICc = zeros(n_models,n_cells);

for c = 1:n_cells
    for m = 1:n_models
        load_model(metas{m,c}.modelpath);
        calc_xxx(1);
        train_name = XXX{end}.training_set{1};
        s = size(XXX{end}.dat.(train_name).stim);
        n = s(1) * s(2);
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
        
        AIC(m,c) = n * log(XXX{end}.score_train_mse) + 2*K;
        AICc(m,c) = n * log(XXX{end}.score_train_mse) + 2*K + (2*K*(K+1))/(n-K-1);
        
    end
end

for c = 1:n_cells
    for m = 1:n_models
%         fprintf(1,'AIC %s = %f\n', modelnames{m,c}, AIC(m,c));
        fprintf(1,'%s %s %f\n', cellids{c} , modelnames{m}, AIC(m,c));
    end
end

end


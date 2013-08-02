function update_batch_files(batch, cellids, modelnames)

global STACK XXX META;

n_models = length(modelnames);
n_cellids = length(cellids);
stacks = cell(n_models, n_cellids);
metas = cell(n_models, n_cellids);
x0s = cell(n_models, n_cellids);

dbopen;
for mi = 1:n_models
    model = modelnames{mi};
    for ci = 1:n_cellids
        cellid = cellids{ci};
        
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(batch) ''];
        sql = [sql ' AND cellid="' cellid '"'];
        sql = [sql ' AND modelname="' model '"'];
        results = mysql(sql);  
                
        fprintf('Updating Model [%d/%d]\n', (mi-1)*n_cellids+ci, n_models*n_cellids);        
        if isempty(results)
            stacks{mi, ci} = {};
            metas{mi, ci} = {};
            x0s{mi, ci} = {};
            fprintf('Skipping because it doesn''t exist.\n');
            continue;
        elseif (length(results) > 1)
            error('Multiple DB Results found: for %s/%s\n', cellid, model);
        end
        
        load_model(char(results(1).modelpath));
        
        calc_xxx(1);
        
        % TODO: REPLACE THIS HARD-CODED STUFF WITH A CALL TO
        % MEASURE_WITH_ALL_METRICS() or something like that, which should
        % also be used when inserting into the DB (after the save)
        
        if isfield(XXX{end}, 'score_test_mse')
            META.perf_est_mse = XXX{end}.score_train_mse;
            META.perf_val_mse = XXX{end}.score_test_mse;            
        end
        
        if isfield(XXX{end}, 'score_test_corr')
            META.perf_est_corr = XXX{end}.score_train_corr;
            META.perf_val_corr = XXX{end}.score_test_corr;
        end
        
        save_model(META.modelpath, STACK, XXX, META);
    end
end
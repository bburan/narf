function update_metrics(batch, cellids, modelnames)

global STACK XXX META MODULES;

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
        
        % -------------------
        
        append_module(MODULES.neural_statistics);        

        % calc_xxx(1);
        
        fns = XXX{1}.training_set;
        avgs = [];
        for ii = 1:length(fns)
            sf = fns{ii};
            avgs(ii) = XXX{end}.dat.(sf).min_dist;
        end        
        META.metric_est_selfdist = sum(avgs); % Mathematically, this isn't true (distances do not always add), but I can't think of a fair way to compare stimuli from different data sets!
        
        fns = XXX{1}.test_set;
        avgs = [];
        for ii = 1:length(fns)
            sf = fns{ii};
            avgs(ii) = XXX{end}.dat.(sf).min_dist;
        end        
        META.metric_val_selfdist = sum(avgs); % Mathematically, this isn't true (distances do not always add), but I can't think of a fair way to compare stimuli from different data sets!
        
        % ---------------------
        save_model(META.modelpath, STACK, XXX, META);
      
    end
end
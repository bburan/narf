function plot_avg_spikes_per_bin_vs_corr(batch, cellids, modelnames)

global STACK XXX META MODULES;

n_models = length(modelnames);
n_cellids = length(cellids);
stacks = cell(n_models, n_cellids);
metas = cell(n_models, n_cellids);
x0s = cell(n_models, n_cellids);

r_test = [];
avg_spikes_per_bin_test = [];
avg_spikes_per_bin_train = [];

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
        r_test(end+1) = results(1).r_test;
        load_model(char(results(1).modelpath));
        calc_xxx(1);

        fns = XXX{1}.test_set;
        v = [];
        %reps = 0;
        for ii = 1:length(fns)
            v = cat(1, v, XXX{end}.dat.(fns{ii}).resp(:));
            %reps = reps + size(XXX{end}.dat.(fns).resp, 3);
        end
        
        avg_spikes_per_bin_test(end+1) = mean(v); 
        
        
        fns = XXX{1}.training_set;
        v = [];
        %reps = 0;
        for ii = 1:length(fns)
            v = cat(1, v, XXX{end}.dat.(fns{ii}).resp(:));
            %reps = reps + size(XXX{end}.dat.(fns).resp, 3);
        end
        
        avg_spikes_per_bin_train(end+1) = mean(v); 
        
        % save_model(META.modelpath, STACK, XXX, META);
        % keyboard;
        
    end
end

figure; 
plot(avg_spikes_per_bin, r_test, 'k.');
xlabel('Spikes per bin');
ylabel('r test');
keyboard;
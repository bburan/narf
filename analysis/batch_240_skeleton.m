function batch_240_skeleton()

batch = 240;
analysis_prefix = 'a240';
baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();   

cells = request_celldb_batch(batch);

mm = {};
mm{1} = module_groups('env100');
mm{2} = module_groups('log2b'); 
mm{3} = module_groups('firn');  
mm{4} = module_groups('npnl', 'senl', 'senl3', 'npfnl', 'npfnl3');  
mm{5} = module_groups('sb', 'fminlsq', 'boost', 'fmin');
mm{6} = module_groups('mse', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6'); 

modulekeys = module_block_combos(mm);

% Enqueue every model
for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        fprintf('Fitting model [%d/%d]\n', jj, length(modulekeys)); 
        fit_single_model(modulekeys{jj}, batch, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
        % TODO: Enqueue into job system instead of doing fit_single_model
        % here to allow work to be distributed everywhere.
        
    end
end

% Generate "top 10" plots for each cellid
% for ii = 1:length(cells)
%     % Plot the top 10 models
%     models = db_get_models(batch, cells{ii}.cellid);
%     compare_models(cellstr(char(models(1:min(length(models), 10)).modelpath)));
%     % TODO: Scatter plots
% end

% TODO: Generate heat map plots (Buggy right now due to caching problem)


% Write a function which makes an SQL query
% Look at db_get_models() for an example
%


function batch_242()

batch = 242;
analysis_prefix = 'a242';
baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();   

cells = request_celldb_batch(batch);

mm = {};

% mm{1} = module_groups('env100');
% mm{2} = module_groups('log2b'); 
% mm{3} = module_groups('firn');  
% mm{4} = module_groups('npnl', 'senl', 'senl3', 'npfnl');  
% mm{5} = module_groups('sb', 'fminlsq', 'boost', 'fmin');
% mm{6} = module_groups('mse', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6'); 

mm{1} = module_groups('env100');
mm{2} = module_groups('log2b'); 
mm{3} = module_groups('firn');  
mm{4} = module_groups('npfnl');  
mm{5} = module_groups('sb', 'fminlsq', 'boost', 'fmin');
mm{6} = module_groups('mse', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6'); 

modulekeys = module_block_combos(mm);

% Enqueue every model
for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        %fprintf('Fitting model [%d/%d]\n', jj, length(modulekeys)); 
        %fit_single_model(modulekeys{jj}, batch, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
        % TODO: Enqueue into job system instead of doing fit_single_model
        % here to allow work to be distributed everywhere.
        enqueue_single_model(modulekeys{jj},  batch, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set, true);
        
    end
end

% Generate "top 10" plots for each cellid
% for ii = 1:length(cells)
%     % Plot the top 10 models
%     models = db_get_models(batch, cells{ii}.cellid);
%     compare_models(cellstr(char(models(max(length(models)-10,1):end).modelpath)));
%     % TODO: Scatter plots
% end

% TODO: Generate heat map plots (Buggy right now due to caching problem)


% Write a function which makes an SQL query
% Look at db_get_models() for an example
%


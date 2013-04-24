function batch_242()

batch = 242;
analysis_prefix = 'a242';
baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();   

cells = request_celldb_batch(batch);

mm = {};

mm{1} = module_groups('env100');
mm{2} = module_groups('log2b'); 
%mm{3} = module_groups('firn');  
%mm{4} = module_groups('npnl', 'senl', 'senl3', 'npfnl');  
%mm{5} = module_groups('sb', 'fminlsq', 'boost', 'fmin');
%mm{6} = module_groups('mse', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6'); 

mm{3} = module_groups('depn');  
mm{4} = module_groups('npfnl');  
mm{5} = module_groups('boost');
mm{6} = module_groups('mses0', 'mses1', 'mses2', 'mses3', 'mses4', 'mses5'); 
force = false;

modulekeys = module_block_combos(mm);

% Enqueue every model
for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        %fprintf('Fitting model [%d/%d]\n', jj, length(modulekeys)); 
        %fit_single_model(modulekeys{jj}, batch, cells{ii}.cellid, ...
        %        cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        % TODO: Enqueue into job system instead of doing fit_single_model
        % here to allow work to be distributed everywhere.
                
        %force = true;
        
        enqueue_single_model(modulekeys{jj},  batch, cells{ii}.cellid, ...
            cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, force);
        
    end
end
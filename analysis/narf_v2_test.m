function narf_v2_test()

baphy_set_path;
narf_set_path; 

global MODULES;

MODULES = scan_directory_for_modules();   

batch = 242;

cells = request_celldb_batch(batch);

mm = {'env100', ...
      'log2b', ...
      'firn', ...
      'npfnl', ...
      {'mse', 'mses5'}, ...
      'boost'};

modulekeys = module_block_combos(mm);

% Enqueue every model
for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        fprintf('Fitting model [%d/%d]\n', jj, length(modulekeys)); 
        
        fit_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
            cells{ii}.training_set, cells{ii}.test_set,  force);
        
%        enqueue_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
%            cells{ii}.training_set, cells{ii}.test_set,  force);
        
    end
end
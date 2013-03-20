function batch_241()
% Enqueues models the following experiment for cells in "batch 241"
% 1. Tries 6 different levels of sparseness penalty
% 2. Tries three different fitters
% 3. Uses the NPNLX module, which lets the NPNL float per-file across the
% training set. The first NPNL is used for the first training & test sets,
% the second NPNL is used for the second training & test sets, and so on.

batch = 241;

baphy_set_path;
narf_set_path; 
global MODULES;
MODULES = scan_directory_for_modules();
cells = request_celldb_batch(batch);

mm = {};
mm{1} = module_groups('env100');
mm{2} = module_groups('log2b');
mm{3} = module_groups('firn');  
mm{4} = module_groups('npnlx');
mm{5} = module_groups('boost', 'sb', 'fminlsq');  
mm{6} = module_groups('mse', 'mses1', 'mses2', 'mses3', 'mses4', 'mses5'); 

modulekeys = module_block_combos(mm); % Combinational fun!
 
for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        % Uncomment next lines to fit models immediately on this machine
        % fprintf('Fitting model [%d/%d]\n', jj, length(modulekeys)); 
        % fit_single_model(modulekeys{jj}, batch, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
        
        % Uncomment next line to enqueue on machine pool
        enqueue_single_model(modulekeys{jj},  batch, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
    end
end

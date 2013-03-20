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
        enqueue_single_model(modulekeys{jj},  batch, ...
            cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
           
    end
end

% Automatically generate unique keys which mask off certain parts of the
% training data set.

% Shift all module groups down one, making space for a new module group in
% the second token position
mm(2:end+1) = mm(1:end);
mm(1) = mm(2);

function ret = query_stephen_about_name(batch, cellid, trainingfile)
    blah = regexp(k, '(\d\d)_', 'tokens');
    ret = blah{1}{1};
end

% Enqueue special models which are trained per-file
for ii = 1:length(cells)
    for jj = 1:length(cells{ii}.training_set)
        % Generate a unique keyname k
        k = ['fm-' query_stephen_about_name(batch,...
                                            cells{ii}.cellid, ...
                                            cells{ii}.training_set{jj})];
        s = [];
        s.(k) = {MODULES.file_masker.mdl(struct('only_indexes', [jj]))}; 
        
        % Insert that new key and module into the 2nd position of mm
        mm{2} = s;
        
        % Then rebuild the combinations and enqueue as before
        modulekeys = module_block_combos(mm);
        
        for kk = 1:length(modulekeys)
            enqueue_single_model(modulekeys{kk},  batch, ...
                cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
        end
    end
end

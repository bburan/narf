function narf_v2_test()

baphy_set_path;
narf_set_path; 

global MODULES;

MODULES = scan_directory_for_modules();   

batch = 242;

mm = {'env100', ...
      'log2b', ...
      {{'firn', 'initrc'}}, ... % {'depn', 'initrc'}, {'add2', 'firn', 'initrc'}, {'firn', 'initrc'}, {'firn', 'init0'}}, ...
      'npfnl', ... % {'npnl', 'npfnl', 'npnlx', 'npfnlx'}, ...
      'mse', ... % 'mse', ... %{'mses0', 'mses1', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6', 'mses7', 'mses8', 'mses9', }, 
      'shboo3'};% {'fmin', 'fminlsq', 'boost', 'fminu', 'qfmin', 'qlsq', 'qboost', 'lsqn', 'genetic', 'anneal', 'sb', 'sp1boost', 'sp2boost', 'sp3boost', 'sp4boost', 'sp5boost'}};

% cells = request_celldb_batch(batch, 'por023b-b1'); % 241
%cells = request_celldb_batch(batch, 'por028d-b1'); % 242
cells = request_celldb_batch(batch);

modulekeys = keyword_combos(mm);

for ii = 1:length(cells)
    for jj = 1:length(modulekeys)       
        fprintf('Fitting cell "%s" [%d/%d] model [%d/%d] %s\n', ...
            cells{ii}.cellid, ii, length(cells),  jj, length(modulekeys), ...
            write_readably(modulekeys{jj})); 
        
        % For testing, use fit_single_model instead of enqueue_single_model
        %fit_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        %  cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        
        enqueue_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, true);
        
    end
end
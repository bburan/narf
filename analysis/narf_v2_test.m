function narf_v2_test()

baphy_set_path;
narf_set_path; 

global MODULES;

MODULES = scan_directory_for_modules();   

% batch = 241;
batch = 242;

mm = {'env100', ...
      'log2b', ...
      {{'add2', 'firn'}, {'depn'}, {'firn'}, {'firn', 'initones'}, {'firn', 'initrc'}}, ...
      {'npnl', 'npfnl'}, ... % {'npnl', 'npfnl', 'npnlx', 'npfnlx'}, ...
      'mse', ...
      'boost'};

%cells = request_celldb_batch(batch, 'por023b-b1');
cells = request_celldb_batch(batch, 'por028d-b1');
%cells = request_celldb_batch(batch);
modulekeys = keyword_combos(mm);

for ii = 1:length(cells)
    for jj = 1:length(modulekeys)       
        fprintf('Fitting cell "%s" [%d/%d] model [%d/%d] %s\n', ...
            cells{ii}.cellid, ii, length(cells),  jj, length(modulekeys), ...
            write_readably(modulekeys{jj})); 
        
        % For testing, use fit_single_model instead of enqueue_single_model
        %fit_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        %    cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        
        %enqueue_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        % cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, true);
        
    end
end
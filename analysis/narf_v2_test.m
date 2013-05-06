function narf_v2_test()

baphy_set_path;
narf_set_path; 

global MODULES;

MODULES = scan_directory_for_modules();   

batch = 242;

mm = {'env100', ...
      'log2b', ...
      {{'firn', 'init0'}, {'firn', 'initrc'}}, ... % {{'add2', 'firn'}, {'depn'}, {'firn'}, {'firn', 'initrc'}}, ...
      'npfnl', ... % {'npnl', 'npfnl', 'npnlx', 'npfnlx'}, ...
      'mse', ...
      {'sp1boost', 'sp2boost', 'sp3boost', 'sp4boost', 'sp5boost'}};% {'fmin', 'fminlsq', 'boost', 'fminu', 'qfmin', 'qlsq', 'qboost', 'lsqn', 'genetic', 'anneal', 'sb'}};

% cells = request_celldb_batch(batch, 'por023b-b1'); % 241
cells = request_celldb_batch(batch, 'por028d-b1'); % 242
%cells = request_celldb_batch(batch);

modulekeys = keyword_combos(mm);

for ii = 1:length(cells)
    for jj = 1:length(modulekeys)       
        fprintf('Fitting cell "%s" [%d/%d] model [%d/%d] %s\n', ...
            cells{ii}.cellid, ii, length(cells),  jj, length(modulekeys), ...
            write_readably(modulekeys{jj})); 
        
        % For testing, use fit_single_model instead of enqueue_single_model
        %fit_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        %   cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        
        enqueue_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, true);
        
    end
end
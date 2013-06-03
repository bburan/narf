function narf_v2_test()

baphy_set_path;
narf_set_path; 

global MODULES;

MODULES = scan_directory_for_modules();   

%  mm = {'env100', ... % {'inferp1','inferp2','inferp3','inferp4','inferp5','inferp6', 'inferp10','inferp20'}, ...
%        'log2b', ...
%        'firn', ... % {'firn', 'firno', 'depn', 'firnv2'}, ...
%        'initrc', ... 
%        'npnl', ... % {'npnl', 'npfnl', 'npnlx', 'npfnlx', 'nonl'}, ...
%        'mse', ... % 'mse', ... %{'mses0', 'mses1', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6', 'mses7', 'mses8', 'mses9', 'err10', 'err15', 'err20'}, 
%        'boostperfile'}; % {'fmin', 'fminlsq', 'boost', 'fminu', 'qfmin', 'qlsq', 'qboost', 'lsqn', 'genetic', 'anneal', 'sb', 'sp1boost', 'sp2boost', 'sp3boost', 'sp4boost', 'sp5boost'}};

mm = {'isi200', 'log2b', {{'firn'}, {'firn', 'npfnl'}}, {'llinv', 'llgam', 'llexp'}, 'boost'}; % , 'log2b', 'firn', 'bic', 'boost'};

batch = 242;
%cells = request_celldb_batch(batch, 'por023b-b1'); % 241
cells = request_celldb_batch(batch, 'por028d-b1'); % 242
%cells = request_celldb_batch(batch);

modulekeys = keyword_combos(mm);

for ii = 1:length(cells)
    for jj = 1:length(modulekeys)       
        fprintf('Fitting cell "%s" [%d/%d] model [%d/%d] %s\n', ...
            cells{ii}.cellid, ii, length(cells),  jj, length(modulekeys), ...
            write_readably(modulekeys{jj})); 
        
        % For testing, use fit_single_model instead of enqueue_single_model
        fit_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
          cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        
        %enqueue_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        % cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, true);
        
    end
end
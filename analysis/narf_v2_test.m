function narf_v2_test()

baphy_set_path;
narf_set_path; 

global MODULES;

MODULES = scan_directory_for_modules();   

batch = 241;
% batch = 242;

% mm = {'env100', ...
%       'log2b', ...
%       {{'firn', {'initrc', 'initzero'}, 'npfnl'}, ...
%        {'depn', 'npnl'}, 'inex'}, ...
%       {'mse', 'mses5'}, ...
%       'boost'};

% initones

mm = {'env100', 'log2b', 'firn', 'npfnlx', 'mse', 'boost'};

cells = request_celldb_batch(batch, 'por023b-b1');
% cells = request_celldb_batch(batch, 'por028d-b1');
modulekeys = keyword_combos(mm);

for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        fprintf('Fitting cell [%d/%d] model [%d/%d]\n', ...
            ii, length(cells),  jj, length(modulekeys)); 
        
        % For testing, use fit_single_model instead of enqueue_single_model
        fit_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
            cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        
        %  enqueue_single_model(batch, cells{ii}.cellid, modulekeys{jj}, ...
        %   cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode);
        
    end
end
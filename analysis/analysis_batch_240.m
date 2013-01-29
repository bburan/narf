% A script to analyze every cell in BAPHY's batch #240.
function analysis_batch_240(instance_number)
% Right now, instance should be 1,2 or 3 so that different matlab instances
% work on different cellids.

addpath('/home/ivar/matlab/baphy');
addpath('/home/ivar/matlab/narf');
baphy_set_path;
narf_set_path; 
global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();

cells = request_celldb_batch(240);

% Define how groups of modules should be combined to make many models
mm = {}; 
mm{1} = module_groups('env100');
mm{2} = module_groups('nocomp', 'log1', 'log2', 'log3', 'log4');
mm{3} = module_groups('fir', 'firb', 'depfir');
mm{4} = module_groups('nonl', 'npnl', 'sig');
mm{5} = module_groups('fminlsq');  % The fit routine is (for now) stored in the correlation module
mm{6} = module_groups('mse');      % The only option right now is to fit using MSE or MSE+smoothness penalty, so I'm storing that in the MSE module

% For every cellid
 for ii = 1:length(cells)   
     if mod(ii + instance_number,3) == 0
        fit_models(mm, cells{ii}.cellid, ...
                        cells{ii}.training_set, ...
                        cells{ii}.test_set);
     end
                
    % Technically fit_models builds a cache already, but when continuing
    % an analysis after debugging or coding, it helps to rebuild the cache
    % with the following:
    %summarize_cellid(cells{ii}.cellid, true); 
    
    % Generate PNGs showing the top 5 models for each cellid.
    %plot_best_models(cells{ii}.cellid, 5, true);
end

% Finally, display and save heat maps of test and training set performance.
% Display only the MM combinations selected above so that we are not 
% overwhelmed with the full plethora of models saved for a particular
% cellid, as may exist because multiple analyses can use the same cellid.
[models, modelnames] = module_combinations(mm);

summary_files = cellfun(@(x) [NARF_SAVED_ANALYSIS_PATH filesep x.cellid '_summary.mat'], ...
                         cells, 'UniformOutput', false);

plot_summary_heatmaps(summary_files);

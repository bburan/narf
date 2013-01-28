% A script to analyze every cell in BAPHY's batch #240.

narf_set_path; 

% Technically these two aren't needed, but if you are editing module code
% then you need to rescan the modules directory after a change.
global MODULES;
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

% Try every model for each in the batch file, and update cache if possible
for ii = 1:length(cells)   
    fit_models(mm, cells{ii}.cellid, ...
                   cells{ii}.training_set, ...
                   cells{ii}.test_set);
    cache_summaries(cells{ii}.cellid, true);
end

% After the analysis has finished, generate PNGs showing the top 5 models
% for each cellid.
for ii = 1:length(cells)
    plot_best_models(cells{ii}.cellid, 5, true);
end

% Finally, display and save heat maps of test and training set performance.
% Display only the MM combinations selected above so that we are not 
% overwhelmed with the full plethora of models saved for a particular
% cellid, as may exist because multiple analyses can use the same cellid.
[models, modelnames] = module_combinations(mm);

analysis_files = cellfun(@(x) [NARF_SAVED_ANALYSIS_PATH filesep x], ...
                         modelnames, 'UniformOutput', false);
for ii = 1:length(cells)
    
                     
                     
plot_saved_analyses('Test Corr', @(c) 100*getfield(c, 'score_test_corr'), analysis_files);
% plot_saved_analyses('Train Corr', @(c) 100*getfield(c, 'score_train_corr'), analysis_files);
% plot_saved_analyses('Fit Time', @(c) getfield(c, 'fit_time'), analysis_files);
% plot_saved_analyses('Exit Code', @(c) getfield(c, 'exit_code'), analysis_files);

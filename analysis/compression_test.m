% A script to analyze a cell file's performance with various pre-FIR
% compression curves. Executes on every cellid in batch #240. 
% At present, each cell takes approximately 30sec x 54 models = 30 minutes
% With 33 cells to test, this should take approximately 18 hours to execute

global MODULES NARF_SAVED_MODELS_PATH NARF_SAVED_ANALYSIS_PATH;

MODULES = scan_directory_for_modules();

% Use request_celldb_batch() to easily use a batch file
cells = request_celldb_batch(240);

% For each cell/train/test combination presented for that batch file
for ii = 1:length(cells)   
    % Test the most likely model combinations
    cellid = cells{ii}.cellid;
    training_set = cells{ii}.training_set;
    test_set = cells{ii}.test_set;
    test_likely_candidates(cellid, training_set, test_set);
end

% Display heat maps of test and training set performance
analysis_files = dir2cell([NARF_SAVED_ANALYSIS_PATH filesep '*.mat'])
analysis_files = cellfun(@(x) [NARF_SAVED_ANALYSIS_PATH filesep x], analysis_files, 'UniformOutput', false);

plot_saved_analyses('Test Corr', @(c) 100*getfield(c, 'score_test_corr'), analysis_files);
plot_saved_analyses('Train Corr', @(c) 100*getfield(c, 'score_train_corr'), analysis_files);
plot_saved_analyses('Fit Time', @(c) getfield(c, 'fit_time'), analysis_files);
plot_saved_analyses('Exit Code', @(c) getfield(c, 'exit_code'), analysis_files);

% After the analysis has finished, you can compare models in detail with:
% for ii = 1:length(cells)
%     cellid = cells{ii}.cellid;
%     training_set = cells{ii}.training_set;
%     test_set = cells{ii}.test_set;
%     test_likely_candidates(cells.cellid, training_set, test_set);
%     filenames = dir2cell([NARF_SAVED_MODELS_PATH filesep cellid filesep '*.mat']);
%     filenames = cellfun(@(f) [NARF_SAVED_MODELS_PATH filesep cellid filesep f], ...
%                             filenames, 'UniformOutput', false);
%     compare_models(filenames);
% end

% Display the top 5 models for a cell


% % % UNCOMMENT AND RUN THIS WHEN DEVELOPING
cellid = 'por025a-b1';
training_set = {'por025a03_p_SPN'}; 
test_set = {'por025a04_p_SPN'};
filenames = dir2cell([NARF_SAVED_MODELS_PATH filesep cellid filesep '*.mat']);
filenames = cellfun(@(f) [NARF_SAVED_MODELS_PATH filesep cellid filesep f], ...
                        filenames, 'UniformOutput', false);
 
% Careful! For some reason matlab often crashes when you compare many
% models together at the same time.
compare_models(filenames(1:5));

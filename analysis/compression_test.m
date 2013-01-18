% A script to analyze a cell file's performance with various pre-FIR
% compression curves. Displays one summary window for each cellid.  

% Use celldb_wrapper to easily use a batch file


% For each cell/train/test combination presented by the celldb wrapper
global NARF_SAVED_MODELS_PATH;

cellid = 'por028d-b1';
training_set = {'por028d03_p_SPN'}; 
test_set = {'por028d04_p_SPN'};

test_likely_candidates(cellid, training_set, test_set);

filenames = dir2cell([NARF_SAVED_MODELS_PATH filesep cellid filesep '*.mat']);
filenames = cellfun(@(f) [NARF_SAVED_MODELS_PATH filesep cellid filesep f], ...
                    filenames, 'UniformOutput', false);

compare_models(filenames);
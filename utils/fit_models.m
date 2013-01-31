function fit_models(mm, cellid, training_set, test_set, skipexisting)
% ARGUMENTS:
%   MM      Cell array of structs whose values are cell arrays of modules
%           and whose fieldnames are textual abbreviations of functionality
%           Example: (TODO)
%               MM{1}.env = {module for baphy that loads envelopes}     
%               MM{2}.log = {log module}
%               MM{2}.sqrt = {sqrt module}
%               MM{3}.sig = {FIR + sigmoidal nonlinearity }
%               MM{3}.exp = {FIR + exponential nonlinearity }
%           For the above MM, 4 models would be fit:
%               log+sig, log+exp, sqrt+sig, sqrt+exp 
%   CELLID
%   TRAINING_SET
%   TEST_SET
%   SKIPEXISTING   If set to false, existing models will be overwritten.
%                  Default is true.

global STACK XXX META MODULES...
    NARF_MODULES_PATH ...
    NARF_SAVED_MODELS_PATH ...
    NARF_SAVED_ANALYSIS_PATH;

if nargin < 5
    skipexisting = true;
end

if isempty(MODULES)
    MODULES = scan_directory_for_modules(NARF_MODULES_PATH);
end

[models, modelnames] = module_combinations(mm);
N_models = length(models);

fprintf('Number of models to be fit %d\n', N_models); 

summary = cell(N_models, 1);

mkdir([NARF_SAVED_MODELS_PATH filesep cellid]);
analysis_file = [NARF_SAVED_ANALYSIS_PATH filesep cellid '_summary.mat'];

% Load existing analysis file, so that incomplete analyses may be continued
if exist(analysis_file, 'file') == 2
    fprintf('Loading existing analysis file.\n');
    summary = getfield(load(analysis_file, 'summary'), 'summary'); 
end

for ii = 1:N_models,
    STACK = {};
    XXX = {};
    XXX{1}.cellid = cellid;
    XXX{1}.training_set = training_set;
    XXX{1}.test_set = test_set;
    
    tic;
    
    META.modelname = modelnames{ii};
    META.modelfile = [cellid '_' modelnames{ii} '_' strcat(training_set{:}) '.mat'];
    META.modelpath = [NARF_SAVED_MODELS_PATH filesep cellid filesep META.modelfile];
    
    fprintf('MODEL [%d/%d]: %s\n', ii, N_models, modelnames{ii});
    
    % If the model savefile exists and we aren't skipping, don't fit
    if skipexisting && exist(META.modelpath, 'file') == 2
        fprintf('Skipping because model file exists.\n');
        continue;
    end
      
    % Append modules from each block one at a time
    model = models{ii};
    for jj = 1:length(model),
        fprintf('Auto-initalizing module [%d/%d]\n', jj, length(model));
        append_module(model{jj}); % Calls auto_init for each
    end
    
    % Fit the model using whatever optimization routine it has
    cormod = find_module(STACK, 'correlation');
    META.exit_code = cormod.fitter();
    META.fitter = func2str(cormod.fitter);
    META.fit_time = toc;
    
    summary{ii} = summarize_model();
    
    save_model(META.modelpath, STACK, XXX, META);
    
    save(analysis_file, 'summary');
end
function success = fit_single_model(modulekeys, batch, cellid, training_set, test_set, filecodes)
% TODO: Documentation
% fit_single_model() is used by the queuing system to train a single model
% as a job which may be run on any machine. 
%
% ARGUMENTS:
%    modulekeys     Cell array of keys to be deciphered by module_groups
%                   and converted into an MM struct.
%    batch
%    cellid
%    training_set
%    test_set
%
% RETURNS: 
%    success        True iff everything was fine

global STACK XXX META MODULES...
    NARF_MODULES_PATH ...
    NARF_SAVED_MODELS_PATH;

if isempty(MODULES)
    MODULES = scan_directory_for_modules(NARF_MODULES_PATH);
end

if ~exist([NARF_SAVED_MODELS_PATH filesep cellid], 'dir')
    mkdir([NARF_SAVED_MODELS_PATH filesep cellid]);
end

if ~exist('filecodes','var'),
    filecodes={};
end

success = false;
STACK = {};
XXX = {};
XXX{1}.cellid = cellid;
XXX{1}.training_set = training_set;
XXX{1}.test_set = test_set;
XXX{1}.filecodes = filecodes;

tic;

% Build the modelname
tmp = cellfun(@(n) sprintf('%s_', n), modulekeys, 'UniformOutput', false);
modelname = strcat(tmp{:});
modelname = modelname(1:end-1); % Remove trailing underscore

% Build the mm struct
mm = {};
for ii = 1:length(modulekeys)
    mm{ii} = module_groups(modulekeys{ii});
end
    
% Build the model 
opts = cellfun(@fieldnames, mm, 'UniformOutput', false);
opt_names = cellfun(@(m) m{1}, opts, 'UniformOutput', false);
blocks = cellfun(@getfield, mm, opt_names, 'UniformOutput', false);
model = {};
for jj = 1:length(blocks),
    for kk = 1:length(blocks{jj})
         model{end+1} = blocks{jj}{kk};
    end
end  

META.modelname = modelname;
META.modelfile = [cellid '_' modelname '_' strcat(training_set{:}) '.mat'];
META.modelpath = [NARF_SAVED_MODELS_PATH filesep cellid filesep META.modelfile];

% Verification of DB and modelfile synchronization
modelfile_exists = (exist(META.modelpath, 'file') == 2);
sql=['SELECT * FROM NarfResults WHERE modelname="' modelname '"'...
    ' AND batch=' num2str(batch) ...
    ' AND cellid="' cellid '"'];
dbopen;
db_results=mysql(sql);

% SVD: Removed checking for existence of model file and
% NarfResults.  Job management now taken care of at a higher
% level. and fit_single_model should just fit no matter what.
if length(db_results) > 1
    error('Multiple DB hits for batch:%d, cellid: %s, modelfile: %s', batch, cellid, modelname);
else
    % There must not be an existing modelfile or DB entry if we reach here. 
    fprintf('Training model, since modelfile not found\n');
    
    % Append modules from each block one at a time
    for jj = 1:length(model),
        fprintf('Auto-initalizing module [%d/%d]\n', jj, length(model));
        append_module(model{jj}); % Calls auto_init for each
    end
    
    % Fit the model using whatever optimization routine it has
    cormod = find_module(STACK, 'correlation');
    META.exit_code = cormod.fitter();
    META.fitter = func2str(cormod.fitter);
    META.fit_time = toc;
    META.batch = batch;
    
    verify_model_polarity(); % invert the model    
    save_model(META.modelpath, STACK, XXX, META);
    db_insert_model();
end

success = true;

end
function success = fit_single_model(batch, cellid, modulekeys, fitterkeys, training_set, test_set, strict_git_logging)
% success = fit_single_model(batch, cellid, modulekeys, fitterkeys, training_set, test_set, strict_git_logging)
%
% Fits a single model described by modul
% 
% fit_single_model() is used by the queuing system to train a single model
% as a job which may be run on any machine. 
% SVD: Removed checking for existence of model file and
% NarfResults.  Job management now taken care of at a higher
% level. and fit_single_model should just fit no matter what.
%
% ARGUMENTS:
%    modulekeys     Cell array of keys to be deciphered by module_groups
%                   and converted into an MM struct.
%    batch          Batch number.
%    cellid         Cell ID.
%    training_set   A cell array of respfiles to fit the model on.
%    test_set       A cell array of respfiles to verify predictive performance.
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

if ~exist('strict_git_logging','var'),
    strict_git_logging = false; % TODO: Change this to true when development stable
end

success = false;
META = [];
STACK = {};
XXX = {};
XXX{1}.cellid = cellid;
XXX{1}.training_set = training_set;
XXX{1}.test_set = test_set;
XXX{1}.filecodes = filecodes;

git = ['git --git-dir=' NARF_PATH '/.git --work-tree=' NARF_PATH ' '];

if strict_git_logging
    % Complain and throw an error if GIT detects outstanding changes.
    [unix_ret_value, ~] = unix([git 'diff-files --quiet']);
    
    if unix_ret_value ~= 0
        fprintf('CMD:%s\n', [git 'diff-files --quiet']);
        error(['\n\n--------------------------------------------------\n' ...
            'ERROR: Unstaged or uncommited changes to NARF detected! \n' ...
            'This is not allowed! We need to store the git commit hash\n' ...
            'to properly mark when a model was fit, and in what exact manner.\n' ...
            'Please commit your changes and try again. Local commits are fine.\n' ...
            '(There is no need to commit your changes to the public repository yet.)\n' ...
            '--------------------------------------------------\n']);
    end
end

cmd = [git 'rev-parse HEAD'];
[~, unix_string] = unix(cmd);
META.git_commit  = regexprep(unix_string, '\n', '');

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
sql=['SELECT * FROM NarfResults WHERE modelname="' modelname '"'...
    ' AND batch=' num2str(batch) ...
    ' AND cellid="' cellid '"'];
dbopen;
db_results=mysql(sql);

if length(db_results) > 1
    error('Multiple DB hits for batch:%d, cellid: %s, modelfile: %s', batch, cellid, modelname);
else
    fprintf('Training model...\n');
    
    % Append modules from each block one at a time
    for jj = 1:length(model),
        fprintf('Auto-initalizing module [%d/%d]\n', jj, length(model));
        append_module(model{jj}); % Calls auto_init for each module
    end
    
    % Fit the model using whatever optimization routine it has
    [cormod, ~] = find_modules(STACK, 'correlation', true);
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
function success = fit_single_model(batch, cellid, modulekeys, training_set, test_set, filecodes, strict_git_logging)
% success = fit_single_model(batch, cellid, modulekeys, training_set, test_set, filecodes, strict_git_logging)
%
% Fits a single model described by modulekeys
% 
% fit_single_model() is used by the queuing system to train a single model
% as a job which may be run on any machine. Doesn't check for the existence
% of a model file or presence in the database, so should not be used for 
% job management.
%
% ARGUMENTS:
%    batch          Batch number.
%    cellid         Cell ID.
%    modulekeys     Cell array of keys to be deciphered by module_groups
%                   and converted into an MM struct.
%    training_set   A cell array of respfiles to fit the model on.
%    test_set       A cell array of respfiles to verify predictive performance.
%    strict_git_logging    When true, throws an error if GIT has any non-
%                          commited files in the repository.
%
% RETURNS: 
%    success        True iff everything was fine

global STACK XXX NARF_PATH META MODULES...
    NARF_MODULES_PATH ...
    NARF_SAVED_MODELS_PATH;

if isempty(MODULES)
    MODULES = scan_directory_for_modules(NARF_MODULES_PATH);
end

if ~exist([NARF_SAVED_MODELS_PATH filesep num2str(batch) filesep cellid], 'dir')
    mkdir([NARF_SAVED_MODELS_PATH filesep num2str(batch) filesep cellid]);
end

if ~exist('strict_git_logging','var'),
    strict_git_logging = false; % TODO: Change this to true when development is stable
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
        fprintf(['\n\n--------------------------------------------------\n' ...
            'ERROR: Unstaged or uncommited changes to NARF detected! \n' ...
            'This is not allowed! We need to store the git commit hash\n' ...
            'to properly mark when a model was fit, and in what exact manner.\n' ...
            'Please commit your changes and try again. Local commits are fine.\n' ...
            '(There is no need to commit your changes to the public repository yet.)\n' ...
            '--------------------------------------------------\n']);
        return;
    end
end

cmd = [git 'rev-parse HEAD'];
[~, unix_string] = unix(cmd);
META.git_commit  = regexprep(unix_string, '\n', '');

% Build the modelname
tmp = cellfun(@(n) sprintf('%s_', n), modulekeys, 'UniformOutput', false);
modelname = strcat(tmp{:});
modelname = modelname(1:end-1); % Remove trailing underscore

META.batch = batch;
META.modelname = modelname;
META.modelfile = [num2str(batch) '_' cellid '_' modelname '.mat'];
META.modelpath = [NARF_SAVED_MODELS_PATH filesep num2str(batch) ...
                  filesep cellid filesep META.modelfile];

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
    
    tic;
    % Build and train the model
    for ii = 1:length(modulekeys)
        run_keyword(modulekeys{ii});
    end
    META.fit_time = toc;    

    verify_model_polarity(); % invert the model    
    save_model(META.modelpath, STACK, XXX, META);
    db_insert_model();
end

success = true;

end
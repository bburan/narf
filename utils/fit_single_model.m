function success = fit_single_model(batch, cellid, keywords_to_exec, training_set, test_set, filecodes, strict_git_logging)
% success = fit_single_model(batch, cellid, keywords_to_exec, training_set, test_set, filecodes, strict_git_logging)
%
% Fits a single model described by keywords_to_exec
% 
% fit_single_model() is used by the queuing system to train a single model
% as a job which may be run on any machine. Doesn't check for the existence
% of a model file or presence in the database, so should not be used for 
% job management.
%
% ARGUMENTS:
%    batch          Batch number.
%    cellid         Cell ID.
%    keywords_to_exec     Cell array of keywords to be executed
%    training_set   A cell array of respfiles to fit the model on.
%    test_set       A cell array of respfiles to verify predictive performance.
%    strict_git_logging    When true, throws an error if GIT has any non-
%                          commited files in the repository.
%
% RETURNS: 
%    success        True iff everything was fine

global STACK XXX NARF_PATH META MODULES...
    NARF_MODULES_PATH NARF_SAVED_MODELS_PATH;

% New default is to scan every time for changes. This prevents stale
% caching problem with memoization.
MODULES = scan_directory_for_modules(NARF_MODULES_PATH);

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
XXX{1}.test_set = {};
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
tmp = cellfun(@(n) sprintf('%s_', n), keywords_to_exec, 'UniformOutput', false);
modelname = strcat(tmp{:});
modelname = modelname(1:end-1); % Remove trailing underscore

META.batch = batch;
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
    % Ivar says: Turning off memoization because disk activity too heavy
    % for the high-resolution data sets that I am doing. March 27, 2014. 
    % Ivar now says: Turning memoization back on. Now it only memoizes
    % functions that took longer than 300 seconds (5 minutes) to run. 
    %memoized_run_keyword = @run_keyword;               % Uncomment for no memoization
    memoized_run_keyword = memoize(@run_keyword);    % Uncomment for memoization
       
    % During development, we actually want to also check that no files have
    % changed before using memoization. It's too hard to figure out if a
    % fitter has changed, or a module has changed, or a util has changed, 
    % or whatever to do this automatically, so we'll be conservative and
    % make sure that absolutely nothing has changed.
    [~, filehash] = unix(['tar -cP ' NARF_PATH ' --exclude .git --to-stdout | md5sum']);
    fprintf('Computing hash for all NARF files: %s\n', filehash)
    
    for ii = 1:length(keywords_to_exec)
        if isempty(findstr(keywords_to_exec{ii},'fit05')),
            % NON-MEMOIZED VERSION:
            run_keyword(keywords_to_exec{ii}, STACK, XXX, META);
        else
            % MEMOIZED VERSION
            fprintf('Calling Keyword: %s\n', keywords_to_exec{ii});
            
            [so, xo, mo] = memoized_run_keyword(keywords_to_exec{ii}, ...
                                                filehash, STACK, XXX, META);
            STACK = so;
            XXX = xo;
            META = mo;
        end
    end
    
    % SVD hacked in from nmse
    % IVAR TODO: Should there be a metric package check here instead?
    mods = find_modules(STACK, 'correlation', true);
    if isempty(mods)
        append_module(MODULES.correlation);    
    end
       
    % Now model specifics occur AFTER all training is done, so that models
    % otherwise similar can use the same memoized files.
    META.modelname = modelname;
    META.modelfile = [num2str(batch) '_' cellid '_' modelname '.mat'];
    META.modelpath = [NARF_SAVED_MODELS_PATH filesep num2str(batch) ...
                  filesep cellid filesep META.modelfile];

    META.fit_time = toc;

    XXX{1}.test_set = test_set;
    XXX = XXX(1);
    calc_xxx(1);
        
    % invert the polarity of the model, possibly screwing it completely up,
    % because there is probably some sort of nonlinearity after the filter
    verify_model_polarity();
    save_model(META.modelpath, STACK, XXX, META);
    db_insert_model();
end

success = true;


end
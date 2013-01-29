function s = summarize_model()
% Returns a structure summarizing the contents of STACK, XXX, and META. 
% This function is often used to build an analysis summary cache file so 
% that the performance and coefficients of many models may be compared
% quickly without having to open each file each time.

global STACK XXX META NARF_PATH;

if length(STACK) > length(XXX)
    recalc_xxx(1);
end

s = [];

% Extract summary of XXX
s.cellid = XXX{1}.cellid;
s.training_set = XXX{1}.training_set;
s.test_set = XXX{1}.test_set;
s.score_train_corr = XXX{end}.score_train_corr;
s.score_test_corr = XXX{end}.score_test_corr;
s.score_train_mse = XXX{end}.score_train_mse;
s.score_test_mse = XXX{end}.score_test_mse;

% Extract summary of STACK
s.n_free_params = length(pack_fittables(STACK));    
s.fits = [];
for ii = 1:length(STACK)
    if ~isfield(STACK{ii}, 'fit_fields')
        continue;
    end
    ff = STACK{ii}.fit_fields;
    for jj = 1:length(ff)
        s.fits.(STACK{ii}.name).(ff{jj}) = STACK{ii}.(ff{jj});
    end
end

% Extract summary of META
s.modelname   = META.modelname;
s.modelfile   = META.modelfile;
s.modelpath   = META.modelpath;
s.fit_time    = META.fit_time;
s.fitter      = META.fitter;
s.exit_code   = META.exit_code;

% Complain and throw an error if GIT detects outstanding changes.
git = ['git --git-dir=' NARF_PATH '/.git ' ...
           '--work-tree=' NARF_PATH ' '];

[unix_ret_value, unix_string] = unix([git 'diff-files --quiet']);

if unix_ret_value ~= 0 
    fprintf('CMD:%s\n', [git 'diff-files --quiet']);
    error(sprintf(['\n\n--------------------------------------------------\n' ...
           'ERROR: Unstaged or uncommited changes to NARF detected! \n' ...
           'This is not allowed! We need to store the git commit hash\n' ...
           'to properly mark when a model was fit, and in what exact manner.\n' ...
           'Please commit your changes and try again. Local commits are fine.\n' ...
           '(There is no need to commit your changes to the public repository yet.)\n' ...
           '--------------------------------------------------\n']));
end

% Otherwise, store the git commit hash
cmd = [git 'rev-parse HEAD'];
[unix_ret_value, unix_string] = unix(cmd);
s.git_commit  = regexprep(unix_string, '\n', '');
function run_keyword(keyword)
% run_keyword(keyword)
% 
% NARF keyword functions are NOT supposed to be on the path, so when you 
% want to run one, you need to use this function. Example:
%
% run_keyword('env100') 
%
global NARF_KEYWORDS_PATH;

warning off MATLAB:dispatcher:nameConflict;
run([NARF_KEYWORDS_PATH filesep keyword]);
warning on MATLAB:dispatcher:nameConflict;
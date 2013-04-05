function run_keyword(keyword)
% run_keyword(keyword)
% 
% NARF keyword functions are NOT supposed to be on the path, so when you 
% want to run one, you need to use this function. Example:
%
% run_keyword('env100') 
%
global NARF_KEYWORDS_PATH;

run([NARF_KEYWORDS_PATH filesep keyword]); 
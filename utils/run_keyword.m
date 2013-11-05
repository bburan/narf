function [stackout, xxxout, metaout] = run_keyword(keyword, ...
     filehash, stack, xxx, meta)
% run_keyword(keyword)
% 
% NARF keyword functions are NOT supposed to be on the path, so when you 
% want to run one, you need to use this function. Example:
%
% Note how this function accepts a bunch of extra arguments. These are
% because now we use a caching system. The hack includes not only the
% keyword, but a MD5 hash of the entire NARF directory so that 
% changes in ANY file can be tracked as if it were a pure function. 
%
% FIXME: It would be nice to keep the old, programmatic interface around,
% so that even if we aren't using memoization, we can still use
% run_keyword.
%
global NARF_KEYWORDS_PATH XXX STACK META;

warning off MATLAB:dispatcher:nameConflict;
run([NARF_KEYWORDS_PATH filesep keyword]);
warning on MATLAB:dispatcher:nameConflict;

stackout = STACK;
xxxout = XXX;
metaout = META;

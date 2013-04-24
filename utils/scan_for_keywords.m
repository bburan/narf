function keywords = scan_for_keywords(keywords_path)
% keywords = scan_for_keywords(keywords_path)
%
% Scans the provided directory if any (or NARF_KEYWORDS_PATH by default) 
% for keyword functions. 
% 
% Keyword functions are NOT supposed to be on the path, so when you want to
% run one, you need to do:
%
% run_keyword('env100') 
%
% Returns a cell array of keywords that were found.

global NARF_KEYWORDS_PATH;

if ~exist('keywords_path', 'var')
    keywords_path = NARF_KEYWORDS_PATH;
end

keywords = dir2cell([keywords_path filesep '*.m']);

for ii = 1:length(keywords)
    tmp = keywords{ii};
    keywords{ii} = tmp(1:end-2); % Cut off .m at end
end
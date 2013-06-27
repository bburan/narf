function tokens = tokenize_string(s)
% tokens = tokenize_string(s)
%
% Converts a string into a cell array of tokens.
%
% ARGUMENTS:
%    s       A string with tokens separated by underscores.
%
% RETURNS:
%    tokens  A cell array of token strings.
%

tokens = regexp(s, '(^.*?)_|(.*?)_|(.*?)$', 'tokens');

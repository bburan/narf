function tokens = tokenize_string(s)
% Converts a modelname into a cell array of tokens
tokens = regexp(s, '(^.*?)_|(.*?)_|(.*?)$', 'tokens');

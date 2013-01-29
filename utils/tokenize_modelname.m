function tokens = tokenize_modelname(modelname)
% Converts a modelname into a cell array of tokens
tokens = regexp(modelname, '(^.*?)_|(.*?)_|(.*?)$', 'tokens');

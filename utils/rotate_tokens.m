function r = rotate_tokens(s, n)
% Rotate a string S so that token number N is first.
% Used in conjunction with sort, you can sort by token easily with this.

tokens = tokenize_modelname(s);

r = strcat(strcat(tokens{n:end}), strcat(tokens{1:n-1}));
r = r{:};
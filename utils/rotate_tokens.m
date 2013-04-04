function r = rotate_tokens(s, n)
% r = rotate_tokens(s, n)
%
% 'Rotates' the underscore-separated tokens in string s so that token
% number n is first.
%
% Used in conjunction with sort, you can sort modelnames or fitternames
% by a specific token index easily with this function.
% 
% ARGUMENTS:
%    s       A modelname string such as 'env100_log2_firn_npnl'
%    n       An integer between 1 and the number of tokens in s.
%
% RETURNS:
%    r       The rotated string with underscores removed.
%
% Example:
%    x = rotate_tokens('env100_log2_firn_npnl', 2)
%    Then x should equal 'log2firnnpnlenv100'. 

tokens = tokenize_string(s);

if length(tokens) < n
    error('Token rotation amount too large for string %s', s);
end

if n < 1
    error('Token rotation amount must be 1 or greater.');
end

if n == 1
    r = strcat(tokens{:});
    r = r{:};
    return;
end

r = strcat(strcat(tokens{n:end}), strcat(tokens{1:n-1}));
r = r{:};
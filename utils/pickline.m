function [l, n_linestyles] = pickline(n)
% [l, n_linestyles] = pickline(n)
%
% Returns a linestyle for a particular index. 
% Useful when 'hold on' is set and you need to plot different lines with
% repeated calls to PLOT() but want them to be the proper linestyles.

linestyle = {'-', '--', '-.', ':'};

n_linestyles = length(linestyle);

l = linestyle{mod(n-1, n_linestyles) + 1};

end
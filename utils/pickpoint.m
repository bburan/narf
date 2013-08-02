function p = pickpoint(n)
% p = pickpoint(n)
%
% Returns a pointstyle for a particular index n.
% Useful when 'hold on' is set and you need to plot different points with
% repeated calls to PLOT() but want them to be the proper point styles.

linestyle = {'.', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'};

n_linestyles = length(linestyle);

p = linestyle{mod(n-1, n_linestyles) + 1};

end
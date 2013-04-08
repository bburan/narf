function do_plot_avg_scatter(stack, xxx, field1, field2)
% A visualization function which plots two fields vs each other
% Sorts and groups into blocks such that there are 100 data points plotted.

% OLD WAY: Concatenate all data points into one big list, across ALL files
%d1 = flatten_field(xxx1{end}.dat, field1);
%d2 = flatten_field(xxx2{end}.dat, field2);

x = xxx{end};

[sf, ~, ~] = get_baphy_plot_controls(stack);
dat = x.dat.(sf);  

if ~isequal(size(dat.(field1)), size(dat.(field2)))
    text(0.35, 0.5, 'Cannot scatter plot a multi-channel input.');
    axis([0, 1, 0 1]);
    return;
end
% Sort and average them by groups of 100
D = [dat.(field1)(:) dat.(field2)(:)]; 
D=D(~isnan(D(:,1))&~isnan(D(:,2)),:);
D = sortrows(D);
D = conv_fn(D, 1, @nanmean, ceil(size(D, 1)/100), 0);

plot(D(:,1), D(:,2), 'k.');
axis tight;

end
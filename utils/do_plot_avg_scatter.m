function do_plot_avg_scatter(stack, xxx, field1, field2)
% A visualization function which plots two fields vs each other
% Sorts and groups into blocks of 100. 

% OLD WAY: Concatenate all data points into one big list, across ALL files
%d1 = flatten_field(xxx1{end}.dat, field1);
%d2 = flatten_field(xxx2{end}.dat, field2);

x = xxx{end};

[sf, stim_idx, unused] = get_baphy_plot_controls(stack);
dat = x.dat.(sf);  

if ~isequal(size(dat.(field1)), size(dat.(field2)))
    text(0.35, 0.5, 'Cannot scatter plot a multi-channel input.');
    axis([0, 1, 0 1]);
    return;
end
% Sort and average them by groups of 100
D = [dat.(field1)(:) dat.(field2)(:)]; 
D = sortrows(D);
D = conv_fn(D, 1, @nanmean, 100, 0);

plot(D(:,1), D(:,2), 'k.');
axis tight;

end
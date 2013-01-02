function do_plot_avg_scatter(xxx1, xxx2, field1, field2)
% A visualization function which plots two fields vs each other
% Sorts and groups into blocks of 100. 

% Concatenate all data points into one big list
d1 = flatten_field(xxx1{end}.dat, field1);
d2 = flatten_field(xxx2{end}.dat, field2);

% Sort and average them by groups of 100
D = [d1 d2]; 
D = sortrows(D);
D = conv_fn(D, 1, @mean, 100, 0);

plot(D(:,1), D(:,2), 'k.');
axis tight;

end
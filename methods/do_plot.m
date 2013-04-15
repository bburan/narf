function c = do_plot(xs, ys, paramset_names, xlab, ylab)
% c = do_plot(xs, ys, paramset_names, xlab, ylab)
%
% TODO

if length(xs) ~= length(ys)
    error('do_plot() needs equal length cell arrays xs, ys');
end

colors = {[0,  0,  1], ...
          [0,  1,  0], ...
          [1,  0,  0], ...
          [0,  1,  1], ...
          [1,  0,  1], ...
          [1,  1,  0], ...
          [0,  0,  0.5], ...
          [0,  0.5,  0], ...
          [0.5,  0,  1], ...
          [0,  0.5,  0.5], ...
          [0.5,  0,  0.5], ...
          [0.5, 0.5,   0]};
      
linestyle = {'-', '--', '-.', ':'};

n_colors = length(colors);
n_linestyles = length(linestyle);

leg = {};
handles = [];
for ii = 1:length(xs)
    x = xs{ii};
    y = ys{ii};
    lw = ceil(ii / n_linestyles); 
    ls = linestyle{mod(ii-1, n_linestyles) + 1};
    for jj = 1:size(ys, 2) 
        c = colors{mod(jj-1, n_colors) + 1};
        % Special case: if there is only one signal, make it black
        if 1 == size(ys,2)
            c = [0, 0, 0];
        end
        h = plot(x, y(:, jj), 'Color', c, 'LineStyle', ls, 'LineWidth', lw);        
        handles(end+1) = h;
        leg{end+1} = [paramset_names{ii} '_ch' num2str(jj)];
    end
end

lh = legend(handles, leg{:}, 'Location', 'NorthWest');
set(lh,'Interpreter','none');

xh = xlabel(xlab, 'Interpreter', 'none');
set(xh, 'Units', 'pixels', 'FontWeight', 'bold');
set(xh,'Position', get(xh,'Position') + [0 20 0]); % Move xlabel up

yh = ylabel(ylab, 'Interpreter', 'none');
set(yh, 'Units', 'pixels', 'FontWeight', 'bold');
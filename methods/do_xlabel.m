function do_xlabel(xlab)

xh = xlabel(xlab, 'Interpreter', 'none');
set(xh, 'Units', 'pixels', 'FontWeight', 'bold');
set(xh,'Position', [50 -5 0]); % Move xlabel up

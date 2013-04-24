function do_plot(xxs, xfield, yfield, sel, xlab, ylab)
% do_plot(xxs, xfield, yfield, sel, xlab, ylab)
% xxs is a cell array of default inputs to a module's .fn function
% xfield (what to use as x axis indexes)
% yfield (what to use as y data)
% sel.stimfile must be a single string
% sel.stim_idx must be a scalar
% sel.chan_idxs can be a scalar, vector, or an empty vector to plot everything
% xlab  x label
% ylab  y axis label

if ~isfield(sel, 'stimfile') || ~isfield(sel, 'stim_idx') || ...
        ~isfield(sel, 'chan_idx')
    error('sel has 3 required fields: stimfile, stim_idx, chan_idx');    
end

hold on;
leg = {};
handles = [];
for ii = 1:length(xxs)
    
    % Skip unless all fields are in existence
    if  ~isfield(xxs{ii}.dat, sel.stimfile) || ...
            ~isfield(xxs{ii}.dat.(sel.stimfile), xfield) || ...
            ~isfield(xxs{ii}.dat.(sel.stimfile), yfield)
        continue;
    end
    
    % Otherwise, plot all the selected channels   
    [ls, n_linestyles] = pickline(ii);
    lw = ceil(ii / n_linestyles); 
    
    chansize = 1:size(xxs{ii}.dat.(sel.stimfile).(yfield), 3);
    
    if isempty(sel.chan_idx)
        sel.chan_idx = chansize;
    end
    
    for qq = 1:length(sel.chan_idx)
        jj = sel.chan_idx(qq);
        c = pickcolor(jj);
        
        % Special case: if there is only one signal possible, make it black
        if 1 == chansize
%            if length(xxs) == 1
                 c = [0, 0, 0];
%             else
%                 c = pickcolor(ii);
%                 ls = pickline(1);
%                 lw = 1;
%             end
        end
        h = plot(xxs{ii}.dat.(sel.stimfile).(xfield)(:), ...
                 squeeze(xxs{ii}.dat.(sel.stimfile).(yfield)(:, sel.stim_idx, jj)), ...
                'Color', c, 'LineStyle', ls, 'LineWidth', lw);        
        handles(end+1) = h;
        leg{end+1} = [yfield '/PS' num2str(ii) '/CH' num2str(jj)];
    end
end

lh = legend(handles, leg{:}, 'Location', 'NorthWest'); 
set(lh,'Interpreter', 'none');

do_xlabel(xlab);
do_ylabel(ylab);

axis tight;

hold off;

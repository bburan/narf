function m = correlation(args)
% A function to linearly correlate two signals. 
% 
% The signals are assumed to be two dimensional matrices at this point.
%  Dimension 1 is time
%  Dimension 2 is the stimuli #
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @correlation;
m.name = 'correlation';
m.fn = @do_correlation;
m.pretty_name = 'Correlation';
m.editable_fields = {'input1', 'input2', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.output = 'score_corr';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_inputs;
m.plot_fns{1}.pretty_name = 'Inputs vs Time';

m.plot_fns{2}.fn = @do_plot_scatter;
m.plot_fns{2}.pretty_name = 'Correlation Scatter Plot';

m.plot_fns{3}.fn = @(stack, xxx) do_plot_avg_scatter(xxx, xxx, stack{end}.input1, stack{end}.input2);
m.plot_fns{3}.pretty_name = 'Smoothed Scatter';

% m.plot_gui_create_fn = @create_chan_selector_gui;

function x = do_correlation(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');

    % Build two very long vectors
    V1 = [];
    V2 = [];
    for sf = x.training_set', sf = sf{1};
        V1 = cat(1, V1, x.dat.(sf).(mdl.input1)(:));
        V2 = cat(1, V2, x.dat.(sf).(mdl.input2)(:));
    end
    
    % Compute the correlation
    R = corrcoef(V1,V2);
    if isnan(R(2,1))
        x.(mdl.output) = R(2,1);
    else
        x.(mdl.output) = R(2,1)^2;
    end
    
end

function do_plot_inputs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    plot(dat.(mdl.time), dat.(mdl.input1)(:, stim_idx), 'b-', ...
         dat.(mdl.time), dat.(mdl.input2)(:, stim_idx), 'g-');
    axis tight;
    legend(mdl.input1, mdl.input2);
    
    % Plot the score in the upper left
    themax = max([max(dat.(mdl.input1)(:, stim_idx)), ...
                  max(dat.(mdl.input2)(:, stim_idx))]);
    text(0, themax , sprintf(' R^2: %f', x.(mdl.output)), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');
     
end

function do_plot_scatter(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    plot(dat.(mdl.input1)(:, stim_idx), ...
         dat.(mdl.input2)(:, stim_idx), 'k.');
    axis tight;
end

end
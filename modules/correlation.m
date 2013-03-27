function m = correlation(args)
% A function to linearly correlate two signals. 
% 
% The signals are assumed to be two dimensional matrices:
%  Dimension 1 is time
%  Dimension 2 is the stimuli #
%
% Returns a function module 'm' which implements the MODULE interface.
% See docs/modules.org for more information.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @correlation;
m.name = 'correlation';
m.fn = @do_correlation;
m.pretty_name = 'Correlation';
m.editable_fields = {'input1', 'input2', 'time', 'fitter', 'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.fitter = @fit_fminlsq;
m.train_score = 'score_train_corr';
m.test_score = 'score_test_corr';
m.output = 'corr_score';  % A score must be minimizeable, use '1/r^2'

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_inputs;
m.plot_fns{1}.pretty_name = 'Inputs vs Time';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_scatter(stack, xxx, stack{end}.input1, stack{end}.input2);
m.plot_fns{2}.pretty_name = 'Correlation Scatter Plot';

m.plot_fns{3}.fn = @(stack, xxx) do_plot_avg_scatter(stack, xxx, stack{end}.input1, stack{end}.input2);
m.plot_fns{3}.pretty_name = 'Smoothed Scatter';

% m.plot_gui_create_fn = @create_chan_selector_gui;

function x = do_correlation(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % -------------
    % Compute the training set correlation, ignoring nans
    V1 = [];
    V2 = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        V1 = cat(1, V1, x.dat.(sf).(mdl.input1)(:));
        V2 = cat(1, V2, x.dat.(sf).(mdl.input2)(:));
    end
    R = corrcoef(excise([V1 V2]));
    if isnan(R)
        x.(mdl.train_score) = NaN;
        x.(mdl.output) = NaN;
    else
        x.(mdl.train_score) = R(2,1);
        x.(mdl.output) = 1 - R(2,1);
    end
    
    %---------------
    % Compute the test set correlation, ignoring nans
    V1 = [];
    V2 = [];
    for ii = 1:length(x.test_set),
        sf = x.test_set{ii};
        V1 = cat(1, V1, x.dat.(sf).(mdl.input1)(:));
        V2 = cat(1, V2, x.dat.(sf).(mdl.input2)(:));
    end
    R = corrcoef(excise([V1 V2]));
    if isnan(R)
        x.(mdl.test_score) = NaN;
    else
        x.(mdl.test_score) = R(2,1);
    end
    
end

function do_plot_inputs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
        
    plot(dat.(mdl.time), dat.(mdl.input1)(:, stim_idx), 'b-', ...
         dat.(mdl.time), dat.(mdl.input2)(:, stim_idx), 'g-');
    
    axis tight;
    legend(mdl.input1, mdl.input2);
    
    % Plot the score in the upper left
    themax = max([max(dat.(mdl.input1)(:, stim_idx)), ...
                  max(dat.(mdl.input2)(:, stim_idx))]);
    text(0, themax , sprintf(' Train r: %f\n Test r : %f', ...
        x.(mdl.train_score), x.(mdl.test_score)), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');

end


end
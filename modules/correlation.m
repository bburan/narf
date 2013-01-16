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
m.editable_fields = {'input1', 'input2', 'time', 'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
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
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');

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
        x.(mdl.train_score) = R;
        x.(mdl.output) = R;
    else
        x.(mdl.train_score) = R(2,1)^2;
        x.(mdl.output) = 1/R(2,1)^2;
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
        x.(mdl.test_score) = R;
    else
        x.(mdl.test_score) = R(2,1)^2;
    end
    
end

function do_plot_inputs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    % removed normalization -- SVD 1/9/13
    %s1 = 1/mean(dat.(mdl.input1)(:, stim_idx));
    %s2 = 1/mean(dat.(mdl.input2)(:, stim_idx));
    s1=1;
    s2=1;
    
    plot(dat.(mdl.time), s1.*dat.(mdl.input1)(:, stim_idx), 'b-', ...
         dat.(mdl.time), s2.*dat.(mdl.input2)(:, stim_idx), 'g-');
    axis tight;
    legend(mdl.input1, mdl.input2);
    
    % Plot the score in the upper left
    themax = max([max(s1.*dat.(mdl.input1)(:, stim_idx)), ...
                  max(s2.*dat.(mdl.input2)(:, stim_idx))]);
    text(0, themax , sprintf(' Train r^2: %f\n Test r^2 : %f', x.(mdl.train_score), x.(mdl.test_score)), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');

end


end
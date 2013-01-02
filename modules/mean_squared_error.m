function m = mean_squared_error(args)
% A function to compute the mean squared error between two signals
% 
% The signals are assumed to be two dimensional matrices at this point.
%  Dimension 1 is time
%  Dimension 2 is the stimuli #
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @mean_squared_error;
m.name = 'mean_squared_error';
m.fn = @do_mean_squared_error;
m.pretty_name = 'Mean Squared Error';
m.editable_fields = {'input1', 'input2', 'time', 'error', 'score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.error  = 'error';
m.score  = 'score_mse';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_inputs_and_mse;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';

function x = do_mean_squared_error(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    score = 0;
    
    % Compute the mean squared error
    for sf = x.training_set', sf = sf{1};
        [T S C] = size(x.dat.(sf).(mdl.input1));
        for s = 1:S,
            error = x.dat.(sf).(mdl.input1)(:,s) - ...
                    x.dat.(sf).(mdl.input2)(:,s);
            x.dat.(sf).(mdl.error)(:, s) = error;
            score = score + mean(error.^2);
        end
    end
    
    x.(mdl.score) = score;
    
end

function do_plot_inputs_and_mse(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    plot(dat.(mdl.time), dat.(mdl.input1)(:, stim_idx), 'b-', ...
         dat.(mdl.time), dat.(mdl.input2)(:, stim_idx), 'g-', ...
         dat.(mdl.time), dat.(mdl.error)(:, stim_idx), 'r-');
    axis tight;
    legend(mdl.input1, mdl.input2, mdl.error);
        
    % Plot the score in the upper left
    themax = max([max(dat.(mdl.input1)(:, stim_idx)), ...
                  max(dat.(mdl.input2)(:, stim_idx)), ...
                  max(dat.(mdl.error)(:, stim_idx))]);
    text(0, themax , sprintf(' MSE: %f', x.(mdl.score)), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');
end

end
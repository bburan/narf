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
m.editable_fields = {'input1', 'input2', 'time', 'error', 'output', ...
                     'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.error  = 'error';
m.train_score  = 'score_train_mse';
m.test_score  = 'score_test_mse';
m.output = 'score_train_mse';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_inputs_and_mse;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';

function x = do_mean_squared_error(mdl, x, stack, xxx)
    % Compute the mean squared error of the training set
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    q = flatten_field(x.dat, x.training_set, mdl.input2); 
    train_score = nanmean((p - q).^2);
    
    % Compute the mean squared error of the test set
    ptest = flatten_field(x.dat, x.test_set, mdl.input1);
    qtest = flatten_field(x.dat, x.test_set, mdl.input2); 
    test_score = nanmean((ptest - qtest).^2);
    
    x.(mdl.train_score) = train_score;
    x.(mdl.test_score) = test_score;   
    x.(mdl.output) = train_score;
end

function do_plot_inputs_and_mse(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    hold on;
    do_plot(xouts, mdls{1}.time, mdls{1}.input1, ...
            sel, '', '');
    do_plot(xouts, mdls{1}.time, mdls{1}.input2, ...
            sel, 'Time [s]', 'Prediction & RespAvg [-]');
    hold off;
end

end

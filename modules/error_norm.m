function m = error_norm(args)
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
m.mdl = @error_norm;
m.name = 'error_norm';
m.fn = @do_error_norm;
m.pretty_name = 'Error Norm';
m.editable_fields = {'input1', 'input2', 'time', 'error', 'output', ...
                     'pnorm', 'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.error  = 'error';
m.pnorm   = 1;  % 1 = L1 Norm, 2 is MSE, etc
m.train_score  = 'score_train_norm';
m.test_score  = 'score_test_norm';
m.output = 'score_train_norm';
m.is_perf_metric = true;

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input1, m.input2, m.time};   % Signal dependencies
m.modifies = {m.error, m.train_score, m.test_score, m.output};  % These signals are modified

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_inputs_and_pnorm;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';

function x = do_error_norm(mdl, x)
    % Compute the mean squared error of the training set
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    q = flatten_field(x.dat, x.training_set, mdl.input2); 
    train_score = abs(p - q);
    train_score = train_score(isfinite(train_score)); % Remove NaNs    
    train_score = sum(train_score .^ (mdl.pnorm)) .^ (1/(mdl.pnorm));
    
    % Compute the mean squared error of the test set
    ptest = flatten_field(x.dat, x.test_set, mdl.input1);
    qtest = flatten_field(x.dat, x.test_set, mdl.input2);     
    test_score = abs(ptest - qtest);
    test_score = test_score(isfinite(test_score)); % Remove NaNs
    test_score = sum(test_score.^(mdl.pnorm)).^(1/mdl.pnorm);
    
    x.(mdl.train_score) = train_score;
    x.(mdl.test_score) = test_score;   
    x.(mdl.output) = train_score;
end

function do_plot_inputs_and_pnorm(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    hold on;
    do_plot(xouts, mdls{1}.time, {mdls{1}.input1, mdls{1}.input2}, ...
            sel, 'Time [s]', 'Prediction & RespAvg [-]');
    hold off;
end

end

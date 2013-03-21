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
                     'smoothness_weight', 'sparseness_weight', ...
                     'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.error  = 'error';
m.smoothness_weight = 0.0;
m.sparseness_weight = 0.0;
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

function x = do_mean_squared_error(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    train_score = 0;
    test_score = 0;
    
    % Compute the mean squared error of the training set
%     for ii = 1:length(x.training_set),
%         sf = x.training_set{ii};
%         [T S C] = size(x.dat.(sf).(mdl.input1));
%         for s = 1:S,
%             error = x.dat.(sf).(mdl.input1)(:,s) - ...
%                     x.dat.(sf).(mdl.input2)(:,s);
%             x.dat.(sf).(mdl.error)(:, s) = error;
%             train_score = train_score + nanmean(error.^2);
%         end
%     end
%     
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    q = flatten_field(x.dat, x.training_set, mdl.input2); 
    train_score = nanmean((p - q).^2);
    
    % Compute the mean squared error of the test set
%     for ii = 1:length(x.test_set),
%         sf = x.test_set{ii};
%         [T S C] = size(x.dat.(sf).(mdl.input1));
%         for s = 1:S,
%             error = x.dat.(sf).(mdl.input1)(:,s) - ...
%                     x.dat.(sf).(mdl.input2)(:,s);
%             x.dat.(sf).(mdl.error)(:, s) = error;
%             test_score = test_score + nanmean(error.^2);
%         end
%     end

    ptest = flatten_field(x.dat, x.test_set, mdl.input1);
    qtest = flatten_field(x.dat, x.test_set, mdl.input2); 

    test_score = nanmean((ptest - qtest).^2);
    
    % Add a penalty related to the non-smoothness of the FIR coefs
    firmod = find_module(stack, 'fir_filter');

    x.(mdl.train_score) = train_score;
    x.(mdl.test_score) = test_score;
    %fprintf('%f\t%f\t%f\n', train_score, mdl.smoothness_weight * diff, ...
    %                         mdl.sparseness_weight * nonpeakiness);
    if ~isfield(mdl, 'sparseness_weight')
        mdl.sparseness_weight = 0;
    end
    
    %fprintf('mse metric: %f\n', sparsity_metric(firmod.coefs(:)));
    %fprintf('MSE: %f\n', train_score + (mdl.sparseness_weight * sparsity_metric(firmod.coefs(:))));
    
    x.sparsity = sparsity_metric(firmod.coefs); % FIXME
    x.(mdl.output) = (train_score) + ...
                     mdl.smoothness_weight * smoothness_metric(firmod.coefs(:)) + ...
                     mdl.sparseness_weight * sparsity_metric(firmod.coefs(:));
end

function do_plot_inputs_and_mse(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    plot(dat.(mdl.time), dat.(mdl.input1)(:, stim_idx), 'b-', ...
         dat.(mdl.time), dat.(mdl.input2)(:, stim_idx), 'g-'); %, ...
         %dat.(mdl.time), dat.(mdl.error)(:, stim_idx), 'r-');
    axis tight;
    legend(mdl.input1, mdl.input2); %, mdl.error);
    
    firmod = find_module(stack, 'fir_filter');
    sparsepenalty = mdl.sparseness_weight * sparsity_metric(firmod.coefs);
    
    % Plot the score in the upper left
    themax = max([max(dat.(mdl.input1)(:, stim_idx)), ...
                  max(dat.(mdl.input2)(:, stim_idx))]); %, ...
                  %max(dat.(mdl.error)(:, stim_idx))]);
    text(0, themax , sprintf(' Train MSE: %f\n Test MSE: %f\n Penalty: %f',...
                             x.(mdl.train_score), x.(mdl.test_score), ...
                             sparsepenalty), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');
end

end

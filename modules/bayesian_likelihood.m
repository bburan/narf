function m = bayesian_likelihood(args)
% For displaying the bayesian likelihood of a noise distribution
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @bayesian_likelihood;
m.name = 'bayesian_likelihood';
m.fn = @do_bayesian_likelihood;
m.pretty_name = 'Bayesian Likelihood';
m.editable_fields = {'stim', 'ISIs', 'time', 'error', 'score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.stim = 'stim';
m.ISIs = 'resp_ISIs';
m.time = 'stim_time';
m.error  = 'error';
m.score  = 'score_mse';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_raw_ISIs;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';

function x = do_mean_squared_error(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    score = 0;
    
    % Compute the mean squared error
    for ii = 1:length(x.training_set) 
        sf = x.training_set{ii};
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
    text(0, themax , sprintf(' BIC: %f', x.(mdl.score)), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');
end

end
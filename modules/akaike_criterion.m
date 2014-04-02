function m = akaike_criterion(args)
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
m.mdl = @akaike_criterion;
m.name = 'akaike_criterion';
m.fn = @akaike_criterion;
m.pretty_name = 'Akaike Information Criterion';
m.editable_fields = {'input1', 'input2', 'time', 'error', 'output', ...
                     'train_score', 'test_score', ...
                     'train_score_norm', 'test_score_norm'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.error  = 'error';
m.train_score  = 'score_akaike';
m.output = 'score_akaike';
m.is_perf_metric = true;

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};

function x = do_akaike_criterion(mdl, x, stack, xxx)
    % Compute the mean squared error of the training set
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    q = flatten_field(x.dat, x.training_set, mdl.input2);
    score_akaike = nanmean((p - q).^2);
    
    x.(mdl.score_akaike) = score_akaike;
    x.(mdl.output) = score_akaike;
end

end

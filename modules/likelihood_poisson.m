function m = likelihood_poisson(args)
% A function to compute the binned likelihood of a poisson process. 
% 
% The signals are assumed to be two dimensional matrices at this point.
%  Dimension 1 is time
%  Dimension 2 is the stimuli #
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @likelihood_poisson;
m.name = 'likelihood_poisson';
m.fn = @do_likelihood_poisson;
m.pretty_name = 'Likelihood Poisson';
m.editable_fields = {'input1', 'input2', 'time', 'output', ...
                     'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.train_score  = 'score_train_nlogl';
m.test_score  = 'score_test_nlogl';
m.output = 'score_train_nlogl';
m.is_perf_metric = true;

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input1, m.input2, m.time};   % Signal dependencies
m.modifies = {m.train_score, m.test_score, m.output};  % These signals are modified

% Optional fields
m.plot_fns = {};
%m.plot_fns{1}.fn = @do_plot_inputs_and_pnorm;
%m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';

function x = do_likelihood_poisson(mdl, x)
    % Compute the mean squared error of the training set
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    r = flatten_field(x.dat, x.training_set, mdl.input2);
    
    train_score = - nansum(r.*log(p) - p);
    
    % Compute the mean squared error of the test set
    ptest = flatten_field(x.dat, x.test_set, mdl.input1);
    rtest = flatten_field(x.dat, x.test_set, mdl.input2);     
    test_score = - nansum(rtest.*log(ptest) - ptest);    
    
    x.(mdl.train_score) = train_score;
    x.(mdl.test_score) = test_score;   
    x.(mdl.output) = train_score;
end

% function do_plot_inputs_and_pnorm(sel, stack, xxx)    
%     [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
%     hold on;
%     do_plot(xouts, mdls{1}.time, {mdls{1}.input1, mdls{1}.input2}, ...
%             sel, 'Time [s]', 'Prediction & RespAvg [-]');
%     hold off;
% end

end

function m = mean_squared_error(args)

global STACK;

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
    'norm_by_se',...
    'train_score', 'test_score', ...
    'train_score_norm', 'test_score_norm'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.time   = 'stim_time';
m.error  = 'error';
m.norm_by_se=0;
m.train_score  = 'score_train_mse';
m.test_score  = 'score_test_mse';
m.train_score_norm  = 'score_train_nmse';
m.test_score_norm  = 'score_test_nmse';
m.output = 'score_train_mse';
m.is_perf_metric = true;

% Overwrite the default module fields with arguments
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input1, m.input2, m.time};   % Signal dependencies
m.modifies = {m.output, m.error};          % These signals are modified

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_inputs_and_mse;
m.plot_fns{1}.fn = @do_plot_inputs_and_mse;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';
m.plot_fns{2}.fn = @do_plot_error_histogram;
m.plot_fns{2}.pretty_name = 'Error Histogram';

    function x = do_mean_squared_error(mdl, x)       
        
        % Compute the mean squared error of the training set
        p = flatten_field(x.dat, x.training_set, mdl.input1);
        q = flatten_field(x.dat, x.training_set, mdl.input2);
        train_score = nanmean((p - q).^2);
        
        % Compute the mean squared error of the test set
        ptest = flatten_field(x.dat, x.test_set, mdl.input1);
        qtest = flatten_field(x.dat, x.test_set, mdl.input2);
        test_score = nanmean((ptest - qtest).^2);
        
        if ~isfield(mdl,'norm_by_se') || ~mdl.norm_by_se,
            train_nmse = train_score / (nanvar(q)+(train_score==0));
            test_nmse = test_score / (nanvar(qtest)+(test_score==0));
        else
            % apply shrinkage filter to nmse
            bincount=10;
            ll=round(linspace(1,length(p)+1,bincount+1));
            llv=round(linspace(1,length(ptest)+1,bincount+1));
            ee=zeros(bincount,1);ve=zeros(bincount,1);
            for bb=1:bincount,
                ee(bb)=nanmean((p(ll(bb):(ll(bb+1)-1))-q(ll(bb):(ll(bb+1)-1))).^2);
                if ~isempty(ptest),
                    ve(bb)=nanmean((ptest(llv(bb):(llv(bb+1)-1))-...
                        qtest(llv(bb):(llv(bb+1)-1))).^2);
                end
            end
            ee=ee./(nanvar(q)+(train_score==0));
            ve=ve./(nanvar(qtest)+(train_score==0));
            me=mean(ee);se=std(ee)./sqrt(bincount);
            train_nmse=shrinkage(me,se,0.5);
            me=mean(ve);se=std(ve)./sqrt(bincount);
            test_nmse=shrinkage(me,se,0.5);
        end
        
        x.(mdl.train_score) = train_score;
        x.(mdl.test_score) = test_score;
        
        if isfield(mdl, 'train_score_norm')
            x.(mdl.train_score_norm) = train_nmse;
            x.(mdl.test_score_norm) = test_nmse;
        end
        x.(mdl.output) = train_score;
    end

    function do_plot_inputs_and_mse(sel, stack, xxx)
        [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end));
        hold on;
        do_plot(xouts, mdls{1}.time, {mdls{1}.input1, mdls{1}.input2}, ...
            sel, 'Time [s]', 'Prediction & RespAvg [-]');
        hold off;
    end

    function do_plot_error_histogram(sel, stack, xxx)
        x = xxx{end};
        mdl = stack{end}{1};
        n_bins = 100;
        
        hist(x.dat.(sel.stimfile).(mdl.input1)(:) - x.dat.(sel.stimfile).(mdl.input2)(:), n_bins);
        
        do_xlabel('STIM Minus RESPAVG [-]');
        do_ylabel('Frequency of Error[-]');
    end

end

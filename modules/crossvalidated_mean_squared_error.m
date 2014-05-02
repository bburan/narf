function m = crossvalidated_mean_squared_error(args)
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
m.mdl = @crossvalidated_mean_squared_error;
m.name = 'crossvalidated_mean_squared_error';
m.fn = @do_crossvalidated_mean_squared_error;
m.pretty_name = 'Cross-validated Mean Squared Error';
m.editable_fields = {'input1', 'input2', 'time', 'error', 'output', ...
    'train_score', 'test_score', ...
    'train_score_norm', 'test_score_norm'};

m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';      % model prediction
m.input2 = 'respavg';   % recorded data
m.time   = 'stim_time';
m.error  = 'error';
m.train_score  = 'score_train_mse';
m.test_score  = 'score_test_mse';
m.train_score_norm  = 'score_train_nmse';
m.test_score_norm  = 'score_test_nmse';
m.output = 'score_train_mse';
m.is_perf_metric = true;

m.crossvalidation_fold = 1;

% Overwrite the default module fields with arguments
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_inputs_and_mse;
m.plot_fns{1}.fn = @do_plot_inputs_and_mse;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';
m.plot_fns{2}.fn = @do_plot_error_histogram;
m.plot_fns{2}.pretty_name = 'Error Histogram';

    function x = do_crossvalidated_mean_squared_error(mdl, x, stack, xxx)
        % Compute the mean squared error of two stratified partitions

        %training_sets = {x.training_set{:}, x.test_set{:}}; % test_set is
        %not available!
        
        training_sets = {x.training_set{:}};
        
        p = flatten_field(x.dat, training_sets, mdl.input1); % model prediction
        q = flatten_field(x.dat, training_sets, mdl.input2); % recorded data
        
        main_partition = p & false;
        c1=1;
        for i=1:length(training_sets),
            duration = size(x.dat.(training_sets{i}).(mdl.input1),1);
            repetition = size(x.dat.(training_sets{i}).(mdl.input1),2);
            main_partition( c1:(c1+duration*ceil(repetition/2)) ) = true;
            c1 = c1+duration*repetition;
        end
        
        other_partition = logical(1-main_partition);
        
        if (mdl.crossvalidation_fold==2),
            tmp = main_partition;
            main_partition = other_partition;
            other_partition = tmp;
        end
        
        train_score = nanmean((p(main_partition) - q(main_partition)).^2);
        train_nmse = train_score / (nanvar(q(main_partition))+(train_score==0));
        
        test_score = nanmean((p(other_partition) - q(other_partition)).^2);
        test_nmse = test_score / (nanvar(q(other_partition))+(test_score==0));
        
        x.(mdl.train_score) = train_score;
        x.(mdl.test_score) = test_score;
        
        if isfield(mdl, 'train_score_norm')
            x.(mdl.train_score_norm) = train_nmse;
            x.(mdl.test_score_norm) = test_nmse;
        end
        x.(mdl.output) = train_score;
        
        % % % WIP: to extend this to K-folds crossvalidation
        %     partitions = [];
        %     c=1;
        %     for i=1:length(x.training_set),
        %         duration = size(x.dat.(x.training_set{i}).(mdl.input1),1);
        %         for j=1:size(x.dat.(x.training_set{i}).(mdl.input1),2),
        %             partitions = [partitions; [c c+duration]];
        %             c = c+duration;
        %         end
        %     end
        %
        %     train_score = 0;
        %     train_nmse = 0;
        %
        %     for i=1:dim(partitions,1),
        %         start_i = partitions(i,1);
        %         end_i = partitions(i,2);
        %         partition_train_score = nanmean((p(start_i:end_i) - q(start_i:end_i)).^2);
        %         partition_train_nmse = partition_train_score / (nanvar(q(start_i:end_i))+(partition_train_score==0));
        %         train_score = train_score + partition_train_score;
        %         train_nmse = train_nmse + partition_train_nmse;
        %     end
        
    end

    function do_plot_inputs_and_mse(sel, stack, xxx)
        [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));
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

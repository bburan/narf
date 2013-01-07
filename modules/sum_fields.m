function m = sum_fields(args)
% Sums multiple fields together. Fields may contain either vectors or
% scalars. 
%
% Useful for collapsing multiple signals into a single signal.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @sum_fields;
m.name = 'sum_fields';
m.fn = @do_sum_fields;
m.pretty_name = 'Sum Fields';
m.editable_fields = {'inputs', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.inputs = {'stim1', 'stim2'}; 
m.time  = 'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_output_vs_time(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_sum_fields(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    % Sum the fields together
    for sf = fieldnames(x.dat)', sf=sf{1};
        x.dat.(sf).(mdl.output) = x.dat.(sf).(mdl.inputs{1});
        for idx = 2:length(mdl.inputs)
            x.dat.(sf).(mdl.output) = x.dat.(sf).(mdl.output) + x.dat.(sf).(mdl.inputs{idx});
        end
    end
end


end
function m = sum_vector_elements(args)
% Sums all the elements in a vector together. Works by default only on the
% third dimension of the input field. 
%
% For collapsing the output channels of a FIR filter into a single signal.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @sum_vector_elements;
m.name = 'sum_vector_elements';
m.fn = @do_sum_vector_elements;
m.pretty_name = 'Sum Vector Elements';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_general_purpose;

% Module fields that are specific to THIS MODULE
m.input  = 'default_stim'; 
m.time   = 'default_stim_time';
m.output = 'default_stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_output_vs_time;
m.plot_fns{1}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_sum_vector_elements(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    for sf = fieldnames(x.dat)', sf=sf{1};
        x.dat.(sf).(mdl.output) = sum(x.dat.(sf).(mdl.input), 3);
    end
end


end
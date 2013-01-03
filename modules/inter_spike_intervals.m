function m = inter_spike_intervals(args)
% Compute the inter-spike intervals (ISIs) of a neural response.
%
% The 1st ISI returned the time between the 1st and 2nd spikes.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @inter_spike_intervals;
m.name = 'inter_spike_intervals';
m.fn = @do_inter_spike_intervals;
m.pretty_name = 'Inter Spike Intervals';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input  = 'resp';
m.time   = 'resp_time';
m.output = 'resp_ISIs';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_raw_ISIs;
m.plot_fns{1}.pretty_name = 'Raw ISIs';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end


function x = do_inter_spike_intervals(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % This analysis will only be valid from the end of the empty space
    % Compute the inter spike intervals referenced from the 'start time'
    % Provide a 'starting time'
    
end

function do_plot_raw_ISIs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    
    hist(x.dat.(sf).(mdl.output_resp)(:), 100);
end

end
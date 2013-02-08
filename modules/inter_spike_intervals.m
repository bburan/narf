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
m.editable_fields = {'input', 'time', 'output', 'output_time'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input  = 'resp';
m.time   = 'resp_time';
m.output = 'resp_ISIs';
m.output_time = 'resp_ISI_times';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_raw_ISIs;
m.plot_fns{1}.pretty_name = 'Plot Raw Spikes';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_inter_spike_intervals(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        in = x.dat.(sf).(mdl.input);
        time =  x.dat.(sf).(mdl.time);
        t_prev = NaN;
        
        [ti, si, ri] = size(in);
        for s = 1:si
            for r = 1:ri
                for t = 1:ti
                    if (t > 1)
                        in(t, s, r)
                    end
                end
            end
        end
    end
end

function do_plot_ISI_spikes(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);
    
    % spikes = dat.(mdl.output)(stim_idx, chan_idx, rep_idx);
    
    % hold on;
    %bar(dat.(mdl.output_time), ones(1,length(spikes)), 0.01,'k');
    % xlabel('Time [s]');
    % ylabel('\lambda(t) [spikes/s]');
    % hold off;
end

end
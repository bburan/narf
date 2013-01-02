function m = normalize_channels(args)
% Normalize RMS power of each channel to be 1.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @normalize_channels;
m.name = 'normalize_channels';
m.fn = @do_normalize_channels;
m.pretty_name = 'Normalize Channels';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time  = 'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_channel_vs_time(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'Channel vs Time';
m.plot_fns{2}.fn = @(stack, xxx) do_plot_channels_as_heatmap(stack, xxx, m.output);
m.plot_fns{2}.pretty_name = 'Channel vs Time (Heatmap)';
m.plot_gui_create_fn = @create_chan_selector_gui;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_normalize_channels(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    chan_sum = [];
    chan_numels = [];
    
    % First work out the power level of every channel
    for sf = fieldnames(x.dat)', sf=sf{1};           
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        
        if isempty(chan_sum)
            chan_sum = zeros(C,1);
            chan_numels = zeros(C,1);
        end
        
        for c = 1:C,
            tmp = x.dat.(sf).(mdl.input)(:,:,c);
            chan_sum(c) = chan_sum(c) + sum(tmp(:));
            chan_numels(c) = chan_sum(c) + numel(tmp);
        end
        
    end
    
    % Compute the RMS levels
    chan_rms = chan_sum ./ chan_numels;
    
    % Now scale every channel by the RMS value
    for sf = fieldnames(x.dat)', sf=sf{1};           
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        
        for c = 1:C,
            x.dat.(sf).(mdl.output)(:,:,c) = 1/chan_rms(c) .* x.dat.(sf).(mdl.input)(:,:,c);
        end
    end
end

end
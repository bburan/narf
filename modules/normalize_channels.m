function m = normalize_channels(args)
% Normalize mean value of each channel to be 1.
% Also, remove DC offset of the channel so that it has zero mean.

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
m.plot_fns{1}.fn = @(stack, xxx) do_plot_all_channels(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'All Normalized Channels';
m.plot_fns{2}.fn = @(stack, xxx) do_plot_channel_vs_time(stack, xxx,  m.time, m.output);
m.plot_fns{2}.pretty_name = 'Single Normalized Channel';
m.plot_fns{3}.fn = @(stack, xxx) do_plot_channels_as_heatmap(stack, xxx, m.output);
m.plot_fns{3}.pretty_name = 'Normalized Channel Heatmap';
m.plot_gui_create_fn = @create_chan_selector_gui;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_normalize_channels(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    % For every channel, remove DC offset and scale by RMS^-1
    for sf = fieldnames(x.dat)', sf=sf{1};        
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        out = zeros(size(x.dat.(sf).(mdl.input)));
        for c = 1:C,
            tmp = x.dat.(sf).(mdl.input)(:,:,c);
            mm = mean(tmp(:));
            rms = sqrt(mean(tmp(:).^2));
            out(:,:,c) = (1/rms) .* (-mm + x.dat.(sf).(mdl.input)(:,:,c));
        end
        x.dat.(sf).(mdl.output) = out;
    end
end

end
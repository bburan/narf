function m = normalize_channels(args)
% Normalize mean value of each channel to be roughly TODO
% Also, remove DC offset of the channel so that it has zero mean.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @normalize_channels;
m.name = 'normalize_channels';
m.fn = @do_normalize_channels;
m.pretty_name = 'Normalize Channels';
m.editable_fields = {'input', 'time', 'output', 'force_positive'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time  = 'stim_time';
m.output = 'stim';
m.force_positive = false;

% Optional fields
%m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_normalized_channels;
m.plot_fns{1}.pretty_name = 'Normalized Channels (All)';
m.plot_fns{2}.fn = @do_plot_single_normalized_channel;
m.plot_fns{2}.pretty_name = 'Normalized Channel (Single)';
%m.plot_fns{3}.fn = @do_plot_normalized_channels_as_heatmap;
%m.plot_fns{3}.pretty_name = 'Normalized Channels (Heatmap)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_normalize_channels(mdl, x, stack, xxx)    
    % Find mean and RMS for all channels across all stimfiles
    tstim=[];
    fns = fieldnames(x.dat);
    for ii = 1:length(fns),
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        tstim = cat(1,tstim,reshape(x.dat.(sf).(mdl.input),T*S,C));
    end
    if mdl.force_positive
        mm = nanmin(tstim);
    else
        mm = nanmean(tstim);
    end
    rms=nanstd(tstim);
    
    % For every channel, remove DC offset and scale by RMS^-1
    fns = fieldnames(x.dat);
    for ii = 1:length(fns),
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        out = zeros(size(x.dat.(sf).(mdl.input)));
        for c = 1:C,
            out(:,:,c) = (1/rms(c)) .* (-mm(c) + x.dat.(sf).(mdl.input)(:,:,c));
        end
        x.dat.(sf).(mdl.output) = out;
    end
end

function do_plot_all_normalized_channels(sel, stack, xxx)       
    [mdls, xins, xouts] = do_calc_paramsets(stack, xxx); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Normalized Signal [-]');
end

function do_plot_single_normalized_channel(sel, stack, xxx)
    [mdls, xins, xouts] = do_calc_paramsets(stack, xxx); 
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Normalized Signal [-]');
end

end
function m = subsample_channels(args)
% Normalize mean value of each channel to be roughly TODO
% Also, remove DC offset of the channel so that it has zero mean.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @subsample_channels;
m.name = 'subsample_channels';
m.fn = @do_subsample_channels;
m.pretty_name = 'Subsample Channels';
m.editable_fields = {'input', 'time', 'output', 'keep_channels'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time  = 'stim_time';
m.output = 'stim';
m.force_positive = false;

% Optional fields
%m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_subsampled_channels;
m.plot_fns{1}.pretty_name = 'Subsampled Chan (All)';
m.plot_fns{2}.fn = @do_plot_single_subsampled_channel;
m.plot_fns{2}.pretty_name = 'Subsampled Chan (Single)';
m.plot_fns{3}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{3}.pretty_name = 'Subsampled Channels (Heatmap)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

function x = do_subsample_channels(mdl, x)    
    % Find mean and RMS for all channels across all stimfiles
    keep_channels=mdl.keep_channels;
    
    % For every channel, remove DC offset and scale by RMS^-1
    fns = fieldnames(x.dat);
    for ii = 1:length(fns),
        sf = fns{ii};
        x.dat.(sf).(mdl.output)=x.dat.(sf).(mdl.input)(:,:,keep_channels);
    
    end        
end

function do_plot_all_subsampled_channels(sel, stack, xxx)       
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Subsampled Signal [-]');
end

function do_plot_single_subsampled_channel(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Subsampled Signal [-]');
end

end

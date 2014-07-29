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
m.plot_fns{1}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{1}.pretty_name = 'Normalized Channels (Heatmap)';
m.plot_fns{2}.fn = @do_plot_single_normalized_channel;
m.plot_fns{2}.pretty_name = 'Normalized Chan (Single)';
m.plot_fns{3}.fn = @do_plot_all_normalized_channels;
m.plot_fns{3}.pretty_name = 'Normalized Chan (All)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

function x = do_normalize_channels(mdl, x)    
    % Find mean and RMS for all channels across all stimfiles
    tstim=[];
    % fns = fieldnames(x.dat); %oops
    fns = x.training_set;
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
               
        out = zeros(T*S,C);
        scale = 1./rms;
        tmp = reshape(x.dat.(sf).(mdl.input), T*S,C);
        for c = 1:C,
            if (rms(c) == 0)
                out(:,c) = tmp(:,c); % div/0 error, just push original data through (doesn't matter anyway)
            else
                out(:,c) = scale(c) .* (-mm(c) + tmp(:,c));
            end
        end
        out = reshape(out,T,S,C);    
        
        if any(isnan(rms))
            % error('divide by zero errors in normalize_channels');
            % We'll just zero everything out so it keeps running and the
            % fitters hopefully realize this is stupid
            x.dat.(sf).(mdl.output) = zeros(size(out));
        else
            x.dat.(sf).(mdl.output) = out;
        end
    end
    
end

function do_plot_all_normalized_channels(sel, stack, xxx)       
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Normalized Signal [-]');
end

function do_plot_single_normalized_channel(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Normalized Signal [-]');
end

end

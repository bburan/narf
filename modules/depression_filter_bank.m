function m = depression_filter_bank(args)
% apply bank of synaptic depression filters to each input channel

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @depression_filter_bank;
m.name = 'depression_filter_bank';
m.fn = @do_depression_filter;
m.pretty_name = 'Depression Filter Bank';
m.editable_fields = {'num_channels', 'strength', 'tau', ...
                     'pedestal', 'input', 'time', 'output' };
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_channels = 2;
m.strength = [0 0.2];  % as fraction of stimulus max magnitude
m.tau = [0 100];  % in ms
m.pedestal = 0;  % as fraction of stimulus max magnitude
m.raw_stim_freq = 100000;
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim_heatmap;
m.plot_fns{1}.pretty_name = 'Filtered Stimulus Heatmap';
m.plot_fns{2}.fn = @do_plot_filtered_stim;
m.plot_fns{2}.pretty_name = 'Filtered Stimulus vs Time';

% Finally, define the 'methods' of this module, as if it were a class
function x = do_depression_filter(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    if ~isfield(mdl,'pedestal'),
        mdl.pedestal=0;
    end
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    %stimmax=x.stimminmax(2,:);
    %keyboard
    
    % calculate global mean level of each channel for scaling dep
    stimmax=[];
    T=0;
    for sf = fieldnames(x.dat)', sf = sf{1};
        if T==0,
            stimmax=nansum(nansum(x.dat.(sf).(mdl.input),1),2);
        else
            stimmax=stimmax+nansum(nansum(x.dat.(sf).(mdl.input),1),2);
        end
        
        T=T+sum(sum(~isnan(x.dat.(sf).(mdl.input)(:,:,1)),1),2);
    end
    stimmax=stimmax./T;
    
    for sf = fieldnames(x.dat)', sf = sf{1};
        [T, S, N] = size(x.dat.(sf).(mdl.input));
        ret = zeros(T, S, N*mdl.num_channels);
        for s = 1:S
            stim_in=squeeze(x.dat.(sf).(mdl.input)(:,s,:))' + ...
                    repmat(mdl.pedestal./stimmax(:),[1,T]);
            depresp=depression_bank(stim_in,(1./stimmax(:))*mdl.strength,...
                                    mdl.tau.*baphy_mod.raw_stim_fs./1000,1)';
            depresp=permute(depresp,[1 3 2]);
            ret(:,s,:) = depresp;
        end
        
        x.dat.(sf).(mdl.output) = ret; 
    end
    
% $$$     % Exotic way to loop over field names using ' and {1}...
% $$$     for sf = fieldnames(x.dat)', sf = sf{1};
% $$$         [T, S, N] = size(x.dat.(sf).(mdl.input));
% $$$         ret = zeros(T, S, N*mdl.num_channels);
% $$$         for s = 1:S
% $$$             stim_in=squeeze(x.dat.(sf).(mdl.input)(:,s,:))';
% $$$             depresp=depression_bank(stim_in,(1./stimmax(:))*mdl.strength,...
% $$$                                     mdl.tau.*baphy_mod.raw_stim_fs./1000,1)';
% $$$             depresp=permute(depresp,[1 3 2]);
% $$$             ret(:,s,:) = depresp;
% $$$         end
% $$$         
% $$$         x.dat.(sf).(mdl.output) = ret; 
% $$$     end
end

% Plot the filter responses
function do_plot_filtered_stim(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    chan_idx = get(baphy_mod.plot_gui.selected_stim_chan_popup, 'Value');
    dat = x.dat.(sf);
    stepsize=baphy_mod.stimulus_channel_count;
    plot(dat.(mdl.time), ...
         squeeze(dat.(mdl.output)(:,stim_idx,chan_idx:stepsize:end)));
    
    axis tight;
end

% Plot heatmap of the filter responses
function do_plot_filtered_stim_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
   % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    chancount=size(dat.(mdl.output),4);
    imagesc(dat.(mdl.time),1:chancount, ...
         squeeze(dat.(mdl.output)(:,stim_idx,:))');
    
    axis tight;
    axis xy
end

end
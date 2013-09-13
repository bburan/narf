
function m = depression_filter_bank(args)
% apply bank of synaptic depression filters to each input channel

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @depression_filter_bank;
m.name = 'depression_filter_bank';
m.fn = @do_depression_filter;
m.pretty_name = 'Depression Filter Bank';
m.editable_fields = {'num_channels', 'strength', 'tau',...
                    'per_channel', 'input', 'time', 'output' };
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_channels = 2;
m.per_channel = 0;
m.strength = [0 0.2];  % as fraction of stimulus max magnitude
m.tau = [0 100];  % in ms
m.tau_norm = 1000;  % 1000 means ms, 10 means 10ths of sec, 1 means sec
m.pedestal = 0;  % as fraction of stimulus max magnitude
m.raw_stim_freq = 100000;
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_depression_cartoon;

m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Filtered Stimuli (All)';
m.plot_fns{2}.fn = @do_plot_single_default_output;
m.plot_fns{2}.pretty_name = 'Filtered Stimuli (Single)';
m.plot_fns{3}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{3}.pretty_name = 'Filtered Stimuli (Heatmap)';
m.plot_fns{4}.fn = @do_depression_cartoon;
m.plot_fns{4}.pretty_name = 'Depression cartoon';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Finally, define the 'methods' of this module, as if it were a class
function x = do_depression_filter(mdl, x, stack, xxx)
    
    if ~isfield(mdl,'pedestal'),
        mdl.pedestal=0;
    end
    if ~isfield(mdl,'tau_norm'),
        mdl.tau_norm=1000;
    end
    [load_mod, ~] = find_modules(stack, 'load_stim_resps_from_baphy',true);
    if isempty(load_mod),
        [load_mod, ~] = find_modules(stack, 'load_stim_resps_wehr',true);
    end
    load_mod = load_mod{1}; % We assume only 1 paramset
    
    % calculate global mean level of each channel for scaling dep
    T=0;
    for sf = fieldnames(x.dat)', sf = sf{1};
        ts=x.dat.(sf).(mdl.input);
        ts=ts.*(ts>0);
        if T==0,
            stimmax=nansum(nansum(ts,1),2);
        else
            stimmax=stimmax+nansum(nansum(ts,1),2);
        end
        T=T+sum(sum(~isnan(ts(:,:,1)),1),2);
    end
    stimmax=stimmax./T;
    [T, S, N] = size(x.dat.(sf).(mdl.input));
    if isfield(mdl,'per_channel') && mdl.per_channel,
        num_channels=length(mdl.tau)./N;
    else
        num_channels=length(mdl.tau);
    end
    
    % Now actually compute depression
    for sf = fieldnames(x.dat)', sf = sf{1};
        [T, S, N] = size(x.dat.(sf).(mdl.input));
        if isfield(mdl,'per_channel') && mdl.per_channel && num_channels<1,
            error('channel count does not work with per_channel==1');
        end
        ret = zeros(T, S, N*num_channels);
        for s = 1:S
            stim_in=squeeze(x.dat.(sf).(mdl.input)(:,s,:))';
            stim_in=stim_in.*(stim_in>0);
            if mdl.pedestal>0,
                thresh = mdl.pedestal./stimmax(:);
                for n=1:N,
                    stim_in(n,stim_in(n,:)<thresh(n))=thresh(n);
                end
                extrabins = round(load_mod.raw_stim_fs*0.5);
                stim_in = cat(2,repmat(thresh,[1 extrabins]),stim_in);
            end
            if isfield(mdl,'per_channel') && mdl.per_channel
                depresp=[];
                for jj=1:N,
                    tresp=depression_bank(...
                        stim_in(jj,:), (1./stimmax(jj))*mdl.strength,...
                        mdl.tau(jj) .* load_mod.raw_stim_fs/mdl.tau_norm, 1);
                    depresp=cat(1,depresp,tresp);
                end
           else
                depresp=depression_bank(...
                    stim_in, (1./stimmax(:))*mdl.strength,...
                    mdl.tau .* load_mod.raw_stim_fs/mdl.tau_norm, 1);
            end
            depresp=permute(depresp',[1 3 2]);
            
            if mdl.pedestal>0,
                depresp = depresp((extrabins+1):end,:,:);
            end
            
            ret(:,s,:) = depresp;
        end
        
        x.dat.(sf).(mdl.output) = ret; 
    end
end

% % Plot the filter responses
function do_depression_cartoon(sel, stack, xxx)
    
    mdl = stack{end}{1};
    fs=stack{1}{1}.raw_stim_fs;
    
    x = struct();
    x.dat.demo.stim=[zeros(fs./2,1);ones(fs,1);zeros(fs,1);ones(fs./5,1);
                 zeros(fs./5,1);ones(fs./5,1);zeros(fs.*1.5,1)];
    xfiltered=do_depression_filter(mdl, x, stack, xxx);
    
    data=[x.dat.demo.stim+1 squeeze(xfiltered.dat.demo.stim)];
    timeaxis=(1:size(data,1))'./fs;
    plot(timeaxis,data);
    
    axis([timeaxis([1 end])' -0.1 2.1]);
    
    legstr=cell(length(mdl.tau)+1,1);
    legstr{1}='Stim';
    for ii=1:length(mdl.tau);
        legstr{ii+1}=sprintf('(u=%.1f,tau=%.3f)',mdl.strength(ii),...
                             mdl.tau(ii)./mdl.tau_norm);
    end
    legend(legstr);
   
end
 
% function do_plot_filtered_stim(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
% 
%     % Find the GUI controls
%     [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
%     chan_idx = get(load_mod.plot_gui.selected_stim_chan_popup, 'Value');
%     dat = x.dat.(sf);
%     stepsize=load_mod.stimulus_channel_count;
%     plot(dat.(mdl.time), ...
%          squeeze(dat.(mdl.output)(:,stim_idx,chan_idx:stepsize:end)));
%     
%     axis tight;
% end
% 
% % Plot heatmap of the filter responses
% function do_plot_filtered_stim_heatmap(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     
%    % Find the GUI controls   
%     [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);
%     chancount=size(dat.(mdl.output),4);
%     imagesc(dat.(mdl.time),1:chancount, ...
%          squeeze(dat.(mdl.output)(:,stim_idx,:))');
%     
%     axis tight;
%     axis xy
% end

end
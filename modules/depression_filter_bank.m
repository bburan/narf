function m = depression_filter_bank(args)
% apply bank of synaptic depression filters to each input channel

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @depression_filter_bank;
m.name = 'depression_filter_bank';
m.fn = @do_depression_filter;
m.pretty_name = 'Depression Filter Bank';
m.editable_fields = {'num_channels', 'strength', 'tau',...
                    'per_channel', 'crosstalk', 'input', 'time', 'output' };
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_channels = 2;
m.per_channel = 0;
m.strength = [0 0.2];  % as fraction of stimulus max magnitude
m.crosstalk = 0;
m.tau = [0 100];  % in ms
m.tau_norm = 1000;  % 1000 means ms, 10 means 10ths of sec, 1 means sec
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
    
    if ~isfield(mdl,'tau_norm'),
        mdl.tau_norm=1000;
    end
    if ~isfield(mdl,'crosstalk'),
        mdl.crosstalk=0;
    end
    [load_mod, ~] = find_modules(stack, 'load_stim_resps_from_baphy',true);
    if isempty(load_mod),
        [load_mod, ~] = find_modules(stack, 'load_stim_resps_wehr',true);
    end
    load_mod = load_mod{1}; % We assume only 1 paramset
    
    % calculate global mean level of each channel for scaling dep
    T=0;
    for sf = fieldnames(x.dat)', 
        sf = sf{1};
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
    stimmax(stimmax==0)=1;
    
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
            
            if isfield(mdl,'per_channel') && mdl.per_channel
                depresp=[];
                if mdl.crosstalk
                    % calc a smoothed "tstim" reflecting crosstalk
                    % that feeds the depression calculator
                    if size(stim_in,1)==2,
                        sfilt=[1-mdl.crosstalk./100 mdl.crosstalk./100;
                               mdl.crosstalk./100 1-mdl.crosstalk./100];
                        ctstim=sfilt*stim_in;
                    elseif size(stim_in,1)>2,
                        sfilt=[mdl.crosstalk./100; ...
                               1-mdl.crosstalk.*2./100; ...
                               mdl.crosstalk./100];
                        ctstim=rconv2(stim_in,sfilt);
                    end
                    for jj=1:N,
                        tresp=depression_bank(...
                            stim_in(jj,:),...
                            (1./stimmax(jj))*mdl.strength(jj)./100,...
                            mdl.tau(jj) .* load_mod.raw_stim_fs/mdl.tau_norm, ...
                            1, ctstim(jj,:));
                        depresp=cat(1,depresp,tresp);
                    end
                else
                    for jj=1:N,
                        tresp=depression_bank(...
                            stim_in(jj,:),...
                            (1./stimmax(jj))*mdl.strength(jj)./100,...
                            mdl.tau(jj).*load_mod.raw_stim_fs/mdl.tau_norm,1);
                        depresp=cat(1,depresp,tresp);
                    end
                end
            else
               depresp=depression_bank(...
                   stim_in, (1./stimmax(:))*mdl.strength./100,...
                   mdl.tau .* load_mod.raw_stim_fs/mdl.tau_norm, 1, ...
                   mdl.crosstalk);
            end
            depresp=permute(depresp',[1 3 2]);
            
            ret(:,s,:) = depresp;
        end
        
        x.dat.(sf).(mdl.output) = ret; 
    end
end

% % Plot the filter responses
function do_depression_cartoon(sel, stack, xxx)
    
    global XXX
    
    mdl = stack{end}{1};
    mdls = stack{end};
    fs=stack{1}{1}.raw_stim_fs;
    tau=[];
    strength=[];
    
    x = struct();
    x.dat.demo.stim=[zeros(fs./2,1);ones(fs,1)./2;zeros(fs,1);
                     ones(round(fs./10),1);zeros(round(fs./10),1);
                     ones(round(fs./10),1);zeros(7.*fs./10,1)];
    if mdl.per_channel,
        x.dat.demo.stim=repmat(x.dat.demo.stim,[1 1 length(mdl.tau)]);
    end
    x.dat.demo.(mdl.input)=x.dat.demo.stim;
    
    if isfield(XXX{1},'filecodes') && ~isempty(XXX{1}.filecodes),
        unique_codes = unique(XXX{1}.filecodes);
    else
        unique_codes={''};
    end
    
    data=[];
    for jj=1:length(mdls),
        xfiltered=do_depression_filter(mdls{jj}, x, stack, xxx);
        xfiltered.dat.demo.stim=xfiltered.dat.demo.(mdl.input);
        
        data=cat(1,data,...
                 [x.dat.demo.stim(:,1)+1 squeeze(xfiltered.dat.demo.stim)]);
        data((end-9):end,:,:)=nan;
    end
    timeaxis=(1:size(data,1))'./fs;
    if size(data,2)>3,
        data=data(:,2:end);
        imagesc(timeaxis,1:size(data,2),data');
        axis xy
        for jj=1:length(mdls),
            %text(1.5+(jj-1)*3.5,2,'stim','VerticalAlign','middle');
            legstr={};
            for ii=1:length(mdls{jj}.tau);
                text(1.55+(jj-1)*3.5,ii,...
                     sprintf('%2d:u=%.1f,{\\tau}=%.0f',...
                             ii,abs(mdls{jj}.strength(ii)),...
                             abs(mdls{jj}.tau(ii)./mdls{jj}.tau_norm.*1000)),...
                     'Color',[1 1 1]);
            end
        end
    else
        plot(timeaxis,data);
        axis([timeaxis([1 end])' -0.1 2.1]);
        for jj=1:length(mdls),
            text(0.5+(jj-1)*3.5,2,[unique_codes{jj} ' stim'],...
                                'VerticalAlign','top');
            legstr={};
            for ii=1:length(mdls{jj}.tau);
                legstr{ii}=sprintf('%2d:u=%.1f,{\\tau}=%.0f',...
                        ii,abs(mdls{jj}.strength(ii)),...
                        abs(mdls{jj}.tau(ii)./mdls{jj}.tau_norm.*1000));
            end
            text(0.5+(jj-1)*3.5,1,legstr,'VerticalAlign','top');
        end
    end
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


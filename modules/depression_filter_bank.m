function m = depression_filter_bank(args)
% apply bank of synaptic depression filters to each input channel

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @depression_filter_bank;
m.name = 'depression_filter_bank';
m.fn = @do_depression_filter;
m.pretty_name = 'Depression Filter Bank';
m.editable_fields = {'num_channels', 'strength', 'tau', 'strength2', 'tau2',...
                    'per_channel', 'offset_in', 'facil_on', 'crosstalk',...
                    'input', 'time', 'output' };
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_channels = 2;
m.per_channel = 0;
m.strength = [0 0.2];  % as fraction of stimulus max magnitude
m.strength2 = [0 0];  % as fraction of stimulus max magnitude
m.tau = [0 100];  % in ms
m.tau2 = [0 0];  % in ms
m.tau_norm = 1000;  % 1000 means ms, 10 means 10ths of sec, 1 means sec
m.offset_in=[0];
m.facil_on = 0;
m.crosstalk = 0;
m.raw_stim_freq = 100000;
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';

% Optional fields
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

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified


% Finally, define the 'methods' of this module, as if it were a class
function x = do_depression_filter(mdl, x)
    
    if ~isfield(mdl,'tau_norm'),
        mdl.tau_norm=1000;
    end
    if ~isfield(mdl,'crosstalk'),
        mdl.crosstalk=0;
    end
    if ~isfield(mdl,'facil_on'),
        mdl.facil_on=0;
    end
    if ~mdl.facil_on,
        mdl.strength=abs(mdl.strength);
    end
    
    % Ivar says: This type of searching can possibly interact with proper
    % operation of the new stack "tree evaluation" rule. For this reason, I
    % broke compatibility with old modules and am removing this. In the
    % future, if the only thing you need to compute is the sampling rate,
    % do that using the appropriate time variable for a signal. 
    %[load_mod, ~] = find_modules(stack, 'load_stim_resps_from_baphy',true);
    %if isempty(load_mod),
    %        [load_mod, ~] = find_modules(stack, 'load_stim_resps_wehr',true);
    %    end
    %    load_mod = load_mod{1}; % We assume only 1 paramset
            
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
        stim_time = x.dat.(sf).(mdl.time);
        raw_stim_fs = 1/(stim_time(2) - stim_time(1)); % A simpler way to get sampling rate
        
        [T, S, N] = size(x.dat.(sf).(mdl.input));
        if isfield(mdl,'per_channel') && mdl.per_channel && num_channels<1,
            error('channel count does not work with per_channel==1');
        end
        ret = zeros(T, S, N*num_channels);
        for s = 1:S
            stim_in=squeeze(x.dat.(sf).(mdl.input)(:,s,:))';
            
            if isfield(mdl,'per_channel') && mdl.per_channel
                % threshold at offset_in level
                if length(mdl.offset_in)==1,
                    stim_in=stim_in-mdl.offset_in;
                    stim_in(stim_in<0)=0;
                else
                    stim_in=stim_in-repmat(mdl.offset_in',[1 T]);
                    stim_in(stim_in<0)=0;
                end

                depresp=[];
                if mdl.crosstalk && size(stim_in,1)>1,
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
                            mdl.tau(jj) .* raw_stim_fs/mdl.tau_norm, ...
                            1, ctstim(jj,:));
                        depresp=cat(1,depresp,tresp);
                    end
                else
                    for jj=1:N,
                        tresp=depression_bank(...
                            stim_in(jj,:),...
                            (1./stimmax(jj))*mdl.strength((jj-1)*num_channels+(1:num_channels))./100,...
                            mdl.tau((jj-1)*num_channels+(1:num_channels)).*raw_stim_fs/mdl.tau_norm,1);
                        depresp=cat(1,depresp,tresp);
                    end
                end
            else
                % threshold at offset_in level
                if length(mdl.offset_in)==1,
                    stim_in=stim_in-mdl.offset_in;
                    stim_in(stim_in<0)=0;
                else
                    stim_in=stim_in-repmat(mdl.offset_in',[1 T]);
                    stim_in(stim_in<0)=0;
                end

                depresp=depression_bank(...
                   stim_in, (1./stimmax(:))*mdl.strength./100,...
                   mdl.tau .* raw_stim_fs/mdl.tau_norm, 1, ...
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
    mdl = stack{end}{1};
    mdls = stack{end};
    
    %fs=stack{1}{1}.raw_stim_fs; % Ivar modified this too   
    for sf = fieldnames(xxx{end}.dat)', sf = sf{1};
        stim_time = xxx{end}.dat.(sf).(mdl.time);
    end
    fs = 1/(stim_time(2) - stim_time(1));
    tau=[];
    strength=[];
    
    x = struct();
    x.dat.demo.(mdl.input)=[zeros(fs./2,1);ones(fs,1)./2;zeros(fs,1);
                     ones(round(fs./10),1);zeros(round(fs./10),1);
                     ones(round(fs./10),1);zeros(7.*fs./10,1)];
    if mdl.per_channel,
        x.dat.demo.(mdl.input)=repmat(x.dat.demo.(mdl.input),[1 1 length(mdl.tau)]);
    end    
    x.dat.demo.(mdl.time)= (1/fs) * (0:size(x.dat.demo.(mdl.input), 1))';
    
    if isfield(xxx{end},'unique_codes') && ~isempty(xxx{end}.unique_codes),
        unique_codes=xxx{end}.unique_codes;
    else
        unique_codes=cell(1,length(mdls));
    end
        
    data=[];
    for jj=1:length(mdls),
        mdls{jj}.offset_in=0;
        x.dat.demo.(mdls{jj}.input)=x.dat.demo.(mdl.input);
        xfiltered=do_depression_filter(mdls{jj}, x);%, stack, xxx);
        xfiltered.dat.demo.stim=xfiltered.dat.demo.(mdls{jj}.input);
        
        data=cat(1,data,...
                 [x.dat.demo.(mdls{jj}.input)(:,1)+1 ...
                  squeeze(xfiltered.dat.demo.(mdls{jj}.input))]);
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
            text(0.5+(jj-1)*3.5,2,[unique_codes{jj} ' ' mdl.input],...
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
    
    do_xlabel('Time [s]');
    do_ylabel('Channel [-]');
end

end


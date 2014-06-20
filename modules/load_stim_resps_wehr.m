function m = load_stim_resps_wehr(args)
% LOAD_STIM_RESPS_WEHR
% A module to load stimulus and response files from Wehr lab-format
% files. Returns a function module which implements the MODULE interface.
% Hacked from load_stim_resps_from_baphy
%

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @load_stim_resps_wehr;
m.name = 'load_stim_resps_wehr';
m.fn = @do_load_wehr;
m.pretty_name = 'Load stim+resp from Wehr Lab data file';
m.editable_fields = {'raw_stim_fs', 'raw_resp_fs', 'include_prestim', ...
                     'stimulus_format','stimulus_channel_count', ...
                     'output_stim', 'output_stim_time', ...
                     'output_resp', 'output_resp_time', 'output_respavg'};
m.isready_pred = @isready_baphy;

% Module fields that are specific to THIS MODULE
m.raw_stim_fs = 100;
m.raw_resp_fs = 100;
m.include_prestim = 1;
m.exclude_target_phase = 1;
m.stimulus_channel_count=0; % 0 should be 'autodetect'
m.stimulus_format = 'envelope';  % Can be 'wav' or 'envelope'
m.output_stim = 'stim';
m.output_stim_time = 'stim_time';
m.output_resp = 'resp';
m.output_resp_time = 'resp_time';
m.output_respavg = 'respavg';
m.is_data_loader = true; % Special marker used by jackknifing routine

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {};   % Signal dependencies
m.modifies = {m.output_stim, m.output_stim_time, m.output_resp, m.output_resp_time ...
                m.output_respavg};          % These signals are modified

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_stim_channels;
m.plot_fns{1}.pretty_name = 'Stim Channels (All)';
m.plot_fns{2}.fn = @do_plot_single_stim_channel;
m.plot_fns{2}.pretty_name = 'Stim Channel (Single)';
% m.plot_fns{3}.fn = @(xx, stck) do_plot_channels_as_heatmap(xx, stck, m.output_stim);
% m.plot_fns{3}.pretty_name = 'Stim Channels (Heatmap)';
m.plot_fns{3}.fn = @do_plot_respavg;
m.plot_fns{3}.pretty_name = 'Response Average';
m.plot_fns{4}.fn = @do_plot_response_raster;
m.plot_fns{4}.pretty_name = 'Response Raster';

% m.plot_fns{6}.fn = @do_plot_stim_log_spectrogram;
% m.plot_fns{6}.pretty_name = 'Stimulus Log Spectrogram';
% m.plot_fns{7}.fn = @do_plot_spectro_and_raster;
% m.plot_fns{7}.pretty_name = 'Spectrogram + Raster';

m.plot_gui_create_fn = @create_gui;

% 
% LOAD INFORMATION ABOUT ALL THE WEHR LAB FILES
%
envpath='/auto/data/daq/wehr/soundfiles/sourcefiles/';
[~,m.file_lookup]=wehr_db;

% ------------------------------------------------------------------------
% Define the 'methods' of this module, as if it were a class

function x = do_load_wehr(mdl, x)
    
    global NARF_DEBUG NARF_DEBUG_FIGURE
    if isempty(NARF_DEBUG),NARF_DEBUG=0;end
    
    % Merge the training and test set names, which may overlap
    files_to_load = unique({x.training_set{:}, x.test_set{:}});
    
    % Create the 'dat' cell array and its entries
    x.dat=struct();
    
    SR=mdl.raw_stim_fs;
    len = length(files_to_load);
    for f_idx = 1:len;
        fbase = files_to_load{f_idx};
        
        % check for special est/val subsets encoded in fbase
        if strcmpi(fbase((end-3):end),'_val'),
            datasubset=2;
            fname=fbase(1:(end-4));
        elseif strcmpi(fbase((end-3):end),'_est'),
            datasubset=1;
            fname=fbase(1:(end-4));
        else
            datasubset=0;
            fname=fbase;
       end
        
        f=m.file_lookup.(fname).file;
        vdim=m.file_lookup.(fname).respfmt;
        if vdim==0,
            vstring='i-clamp';
        elseif vdim==1,
            vstring='v-clamp_E';
        elseif vdim==2,
            vstring='v-clamp_I';
        end
        
        if strcmp(computer,'PCWIN'),
            ppdir=tempdir;
        elseif exist('/auto/data/tmp/tstim/','dir')
            ppdir=['/auto/data/tmp/tstim/'];
        else
            ppdir=['/tmp/'];
        end
        cachefilename=sprintf('%sload_stim_resps_wehr_%s_%s_fs%d.mat',...
                              ppdir,fname,vstring,SR);
        if ~NARF_DEBUG && exist(cachefilename,'file'),
            fprintf('loading cached file %s\n',cachefilename);
            load(cachefilename);
        else
            load(f);
            
            pp=fileparts(f);
            envfile=strrep(out.epochfilenames{1},'\',filesep);
            envfile=strrep(basename(envfile),'sourcefile','envfile');
            e=load([envpath envfile]);
            
            SRint=e.EnvSamplingRate;
            secPerSegment=size(e.EnvSet,1)./length(e.sequence)./SRint;
            
            s0=zeros(SRint,size(e.EnvSet,2));
            
            % stim/resp data will be loaded into s/r matrices
            % SR Hz sampling rate for both
            s=[];r=[];onsettimes=[];
            setcount=size(out.M1,1);
            for ii=1:setcount,
                
                % average response across repetitions
                if size(out.M1,4)>1,
                    % v-clamp data
                    % dim 2 index specifies which v-vclamp mode
                    %ri=squeeze(mean(out.M1(ii,vdim,1:out.nreps(ii,vdim),:),3));
                    ri=squeeze(out.mM1(ii,vdim,:));
                    if vdim==1,
                        % flip sign for E current so that up="more"
                        ri=-ri;
                    end
                    % load low-res stimulus just to make sure it's aligned with
                    % the envelope
                    dsi=squeeze(mean(out.M1stim(ii,1,1:out.nreps(ii),:),3));
                    sentencesperoffset=size(out.sequences,4);
                else
                    % i-clamp data
                    ri=squeeze(mean(out.M1(ii,:,:),2));
                    %ri=squeeze(out.mM1(ii,:))';
                    dsi=squeeze(mean(out.M1stim(ii,1:out.nreps(ii),:),2));
                    sentencesperoffset=size(out.sequences,3);
                end
                ri=resample(ri,SRint,out.samprate);
                dsi=resample(dsi,SRint,out.samprate);
                
                % number of SRint hz samples per segment (segment==sentence?)
                offset=(secPerSegment*SRint)*sentencesperoffset;
                %ii
                %keyboard
                if ii<size(out.M1,1) || offset.*ii<size(e.EnvSet,1),
                    si=[s0;e.EnvSet(offset.*(ii-1)+(1:offset),:);s0];
                else
                    si=[s0;e.EnvSet((offset*(ii-1)):end,:);s0];
                end
                
                if SR<SRint,
                    ri=resample(ri,SR,SRint);
                    dsi=resample(dsi,SR,SRint);
                    si0=si;
                    si=[];
                    for jj=1:size(si0,2);
                        tsi=resample(si0(:,jj),SR,SRint);
                        si=[si tsi];
                    end
                end
                
                % trim possible artifact, last 0.75 sec of r
                ri=ri(1:(end-round(0.75*SR)),:);
                
                % make sure stim and resp are the same length
                if length(ri)>length(si),
                    ri=ri(1:length(si),:);
                    dsi=dsi(1:length(si),:);
                else
                    si=si(1:length(ri),:);
                    dsi=dsi(1:length(ri),:);
                end
                
                %remove possible artifacts from 0.1 sec onset
                ri=ri(round(0.1*SR)+1:end);
                si=si(round(0.1*SR)+1:end,:);
                dsi=dsi(round(0.1*SR)+1:end);
                
                if 1,
                    % remove trend by subtracting min-filtered version
                    % of signal
                    N=round(SR./2);
                    sy=gsmooth(ri(:),N./2);
                    y=zeros(size(sy(:)));
                    for ii=1:length(ri(:));
                        xx=max(1,ii-N);
                        yy=min(length(sy(:)),ii+N);
                        y(ii)=min(sy(xx:yy));
                    end
                    y=gsmooth(y,N./2);
                    
                    if NARF_DEBUG,
                        if isempty(NARF_DEBUG_FIGURE),
                            NARF_DEBUG_FIGURE=figure;
                        end
                        sfigure(NARF_DEBUG_FIGURE);
                        clf
                        subplot(2,1,1);
                        plot([ri(:) sy ri(:)-y y]);
                        axis tight
                        subplot(2,1,2);
                        ty=ri(:)-y;
                        ts=nanmean(si,2);
                        plot([ts./max(ts) ty./max(ty)]);
                        axis tight
                        pause(0.1);
                        %keyboard
                    end
                    ri=ri(:)-y;
                    ri=sqrt(abs(ri)).*sign(ri);
                    ri=ri;
                elseif 0,
                    % remove linear trend
                    ri=detrend(ri);
                else
                    % remove trend by substracting 5th order polynomial fit
                    order = 5;
                    %p = polyfit((1:numel(ri))', ri(:), order);
                    p = polyfit((1:numel(ri))', gsmooth(ri(:),50), order);
                    %plot([ri(:)  gsmooth(ri(:),100) polyval(p, (1:numel(ri))')]);
                    
                    if NARF_DEBUG,
                        if isempty(NARF_DEBUG_FIGURE),
                            NARF_DEBUG_FIGURE=figure;
                        end
                        sfigure(NARF_DEBUG_FIGURE);
                        clf
                        plot([ri(:) polyval(p, (1:numel(ri))') ...
                              ri(:) - polyval(p, (1:numel(ri))')]);
                        keyboard;
                    end
                    
                    ri = ri(:) - polyval(p, (1:numel(ri))');
                end
                
                onsettimesi=(0:floor(length(si)./(secPerSegment*SR)-1))*...
                    (secPerSegment.*SR)+1+round(length(s0)./SRint*SR)-10-f_idx;
                onsettimes=[onsettimes onsettimesi+length(r(:))];
                
                if length(si)<size(s,1),
                    sdiff=size(s,1)-length(si);
                    si=cat(1,si,ones(sdiff,size(si,2)).*nan);
                    ri=cat(1,ri,ones(sdiff,size(ri,2)).*nan);
                end
                si=permute(si,[1 3 2]);
                
                r=[r ri];
                s=[s si];
            end
            fprintf('Saving cache file %s\n',cachefilename);
            esequence=e.sequence;
            save(cachefilename,'r','s','onsettimes','secPerSegment',...
                 'esequence');
        end
        
        if ~mdl.include_prestim,
            % nan-out one-second pre-silence and first second of stim 
            % (ie, get rid of onset transients!)
            r(1:(onsettimes(1)+1*SR),:)=nan;
        end
        
        resp=r;
        stim=s;
        % TODO: Right now the envelope is not necessarily
        % positive semidefinite, and when run through a square root
        % compressor results in complex numbers being developed. This
        % should be fixed on the baphy side, but for now let's just
        % add a workaround here.
        stim(stim<0)=0;
        
        % calculate min and max of all files in the training
        % set... used by depression_filter_bank
        stimminmax=[squeeze(nanmin(nanmin(stim,[],1),[],2))';
                    squeeze(nanmax(nanmax(stim,[],1),[],2))'];
        if f_idx==1,
            x.stimminmax=stimminmax;
        elseif f_idx<=length(x.training_set),
            x.stimminmax=[nanmin(x.stimminmax(1,:),stimminmax(1,:));
                          nanmax(x.stimminmax(2,:),stimminmax(2,:))];
        end
        
        % SVD 2013-03-08 - if specified, pull out either estimation (fit) or
        % validation (test) subset of the data
        if datasubset,
            % break out estimation and validation sets
            % maybe not a hack?
            slen=secPerSegment*SR;  
            
            ffval=find(esequence(1:length(onsettimes))==1);
            ffest=find(esequence(1:length(onsettimes))~=1);
            
            FitRange=[];
            for ffe=ffest,
                FitRange=cat(1,FitRange,onsettimes(ffe)+(1:slen)');
            end
            ValRange=[];
            for ffv=ffval,
                ValRange=cat(1,ValRange,onsettimes(ffv)+(1:slen)');
            end
            
            if datasubset==2
                resp(FitRange)=nan;
             else
                resp(ValRange)=nan;
            end
            
            % nan out the pre-stim silence:
            resp(1:(onsettimes(1)-1),:)=nan;
            % nan out the post-stim silence:
            resp((onsettimes(end)+slen+1):end)=nan;
            resp((onsettimes(max(find(onsettimes<size(resp,1))))+slen+1):end,:)=nan;
        end
        
        x.dat.(fbase).(mdl.output_stim) = stim;
        x.dat.(fbase).(mdl.output_resp) = resp;
        x.dat.(fbase).(mdl.output_respavg) = ...
            squeeze(nanmean(x.dat.(fbase).(mdl.output_resp),3));
        
        % Scale respavg so it is a spike rate in Hz
        x.dat.(fbase).(mdl.output_respavg) = (mdl.raw_resp_fs) .* x.dat.(fbase).(mdl.output_respavg);
        
        % Create time signals for later convenience
        [s1 s2 s3] = size(x.dat.(fbase).(mdl.output_stim));
        [r1 r2 r3] = size(x.dat.(fbase).(mdl.output_resp));
        [a1 a2]    = size(x.dat.(fbase).(mdl.output_respavg));
        
        % TODO: Check stim, resp, raw_time signal sizes match.       
        x.dat.(fbase).resptype = vstring;
        x.dat.(fbase).(mdl.output_stim_time) = (1/mdl.raw_stim_fs).*[1:s1]';
        x.dat.(fbase).(mdl.output_resp_time) = (1/mdl.raw_resp_fs).*[1:r1]';
    end
    
    disp('load_stim_resps_wehr done');
end

% ------------------------------------------------------------------------
% Helper functions

    function ylab = what_is_ylabel(mdls)
        if strcmp(mdls{1}.stimulus_format, 'envelope')
            ylab = 'Envelope Magnitude [?]';
        elseif strcmp(mdls{1}.stimulus_format, 'wav')
            ylab = 'Volume [?]';
        else
            ylab = 'Unknown';
        end
    end

% ------------------------------------------------------------------------
% Plot functions

function do_plot_all_stim_channels(sel, stack, xxx)       
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end));
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);        
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.output_stim_time, mdls{1}.output_stim, ...
            sel, 'Time [s]', what_is_ylabel(mdls));
end

function do_plot_single_stim_channel(sel, stack, xxx)   
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);
    do_plot(xouts, mdls{1}.output_stim_time, mdls{1}.output_stim, ...
            sel, 'Time [s]', what_is_ylabel(mdls)); 
end

function do_plot_respavg(sel, stack, xxx)
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);
    sel.chan_idx = 1;    
    do_plot(xouts, mdls{1}.output_resp_time, mdls{1}.output_respavg, ...
            sel, 'Time [s]', 'Spike Rate Average [Hz]'); 
end

function do_plot_response_raster(sel, stack, xxx)    
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);
    mdl = mdls{1};
    dat = xouts{1}.dat.(sel.stimfile);
    [T, S, R] = size(dat.(mdl.output_resp));
    hold on;
    for r = 1:R
        [xs,ys] = find(dat.(mdl.output_resp)(:, sel.stim_idx, r) > 0);
        plot(dat.(mdl.output_resp_time)(xs), r*ys, 'k.');
    end
    do_xlabel('Time [s]');
    do_ylabel('Trial #');
    axis([0 dat.(mdl.output_resp_time)(end) 0 R+1]);
    hold off;
end

% function do_plot_respavg_as_spikes(stack, xxx)
%     mdl = stack{end};    
%     x = xxx{end};
%     
%     % Read the GUI to find out the selected stim files
%     sf = popup2str(mdl.plot_gui.selected_stimfile_popup);
%     stim = popup2num(mdl.plot_gui.selected_stim_idx_popup);
%     dat = x.dat.(sf);
%     
%     [xs,ys] = find(dat.respavg(idx, :) > 0);
%     bar(dat.(mdl.output_resp_time)(ys), dat.respavg(idx,ys), 0.01, 'k-');
%     axis([0 dat.(mdl.output_resp_time)(end) 0 max(dat.respavg(:, stim))]);
% 
% end

% function do_plot_stim_log_spectrogram(stack, xxx)
%     mdl = stack{end};    
%     x = xxx{end};
%     
%     if strcmp(mdl.stimulus_format, 'envelope')
%         text(0.35, 0.5, 'Cannot visualize envelope as spectrogram');
%         axis([0, 1, 0 1]);
%         return;
%     end
%     
%     % Read the GUI to find out the selected stim files
%     sf = popup2str(mdl.plot_gui.selected_stimfile_popup);
%     stim = popup2num(mdl.plot_gui.selected_stim_idx_popup);
%     chan = popup2num(mdl.plot_gui.selected_stim_chan_popup);
%     
%     dat = x.dat.(sf);
%    
%     % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
%     logfsgram(dat.(mdl.output_stim)(:, stim, chan), ...
%               4048, mdl.raw_stim_fs, [], [], 500, 12); 
%     caxis([-20,40]);  % TODO: use a 'smarter' caxis here
%     axis tight;
% end
% 
% function do_plot_spectro_and_raster(stack, xxx)
%     mdl = stack{end};    
%     x = xxx{end};
%     
%     if strcmp(mdl.stimulus_format, 'envelope')
%         text(0.35, 0.5, 'Cannot visualize envelope as spectrogram');
%         axis([0, 1, 0 1]);
%         return;
%     end
%     
%     % Read the GUI to find out the selected stim files
%     sf = popup2str(mdl.plot_gui.selected_stimfile_popup);
%     stim = popup2num(mdl.plot_gui.selected_stim_idx_popup);
%     dat = x.dat.(sf);
%     
%     hold on;
%     % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
%     logfsgram(dat.(mdl.output_stim)(:,stim)', 4048, mdl.raw_stim_fs, [], [], 500, 12); 
%     caxis([-20,40]);  % TODO: use a 'smarter' caxis here
%     h = get(gca, 'YLim');
%     d = h(2) - h(1);
%     axis tight;    
%     [T, S, R] = size(dat.(mdl.output_resp));
%     hold on;
%     for r = 1:R
%         [xs,ys] = find(dat.(mdl.output_resp)(:, stim, r) > 0);
%         plot(dat.(mdl.output_resp_time)(xs), ...
%              h(1) + (r/R)*d*(d/(d+d/R))*dat.(mdl.output_resp)(xs,stim,r), 'k.');
%     end    
%     hold off;
% end


function hs = create_gui(parent_handle, stack, xxx)
    pos = get(parent_handle, 'Position');
    w = pos(3) - 5;
    h = pos(4) - 5;
    hs = [];
    
    mdl = stack{end}{1};
    mod_idx = length(stack);
    x = xxx{end};
    
    % Create a popup which selects
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'Stimfile:', ...
        'Units', 'pixels', 'Position', [5 (h-25) 100 25]);
    hs.selected_stimfile_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [5 (h-40) w-5 25], ...
        'Callback', @(a,b,c) selected_stimfile_popup_callback());
    
    % Create a stimfile selector
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left', 'String', 'Idx:', ...
        'Units', 'pixels', 'Position', [5 (h-70) 50 25]);
    hs.selected_stim_idx_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'NONE', ...
        'Units', 'pixels', 'Position', [45 (h-70) w-50 25], ...
        'Callback', @(a,b,c) selected_stim_idx_popup_callback());

    % Create a channel selector
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'Chan:', ...
        'Units', 'pixels', 'Position', [5 (h-95) 50 25]);
    hs.selected_stim_chan_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [45 (h-95) w-50 25], ...
        'Callback', @(a,b,c) selected_stim_chan_popup_callback());
    
    hs.textbox = uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', '', ...
        'Units', 'pixels', 'Position', [5 (h-125) w-5 h-100]);
    
    % Two functions to populate the two popup menus
    function update_selected_stimfile_popup()
        fns = fieldnames(x.dat);
        set(hs.selected_stimfile_popup, 'String', char(fns));
        
    end
    
    function update_set_textbox()
        % Print "TEST SET" or "TRAINING SET" in plot GUI as appropriate
        sf = popup2str(hs.selected_stimfile_popup);
        is_test = any(strcmp(sf, x.test_set));
        is_training = any(strcmp(sf, x.training_set));
        
        if is_training & is_test
            str = 'Est&Val Set';
        elseif is_training
            str = 'Estimation Set';
        elseif is_test
            str = 'Validation Set';
        else
            str = 'Not in a set!?'
        end
        
        % Also append the filecode info
        if isfield(x, 'filecodes') && ~isempty(x.filecodes)
            if ~isempty(x.test_set),
                fc = x.filecodes(or(strcmp(sf, x.training_set), ...
                                    strcmp(sf, x.test_set)));
            else
                fc = x.filecodes(strcmp(sf, x.training_set));
            end
            str = [str ' [' sprintf('%s', fc{1}) ']'];
        end
        
        set(hs.textbox, 'String', str); 
    end
    
    function update_selected_stim_idx_popup()
        sf = popup2str(hs.selected_stimfile_popup);
        
        if isfield(x.dat, sf)
            [d1, d2, d3] = size(x.dat.(sf).(mdl.output_stim));
            d = {};
            for i = 1:d2
                d{i} = sprintf('%d',i);
            end
            set(hs.selected_stim_idx_popup, 'String', char(d));
            set(hs.selected_stim_idx_popup, 'Value', 1);
        else
            error('Selected stimulus file not found: %s', sf);
        end
    end
    
    function update_selected_stim_chan_popup()
        sf = popup2str(hs.selected_stimfile_popup);
        
        if isfield(x.dat, sf)
            [d1, d2, d3] = size(x.dat.(sf).(mdl.output_stim));
            d = {};
            for i = 1:d3
                d{i} = sprintf('%d',i);
            end
            set(hs.selected_stim_chan_popup, 'String', char(d));
            set(hs.selected_stim_chan_popup, 'Value', 1);
        else
            error('Selected stimulus file not found: %s', sf);
        end
    end
    
    % Call the two update functions once to build the lists
    update_selected_stimfile_popup();
    update_selected_stim_idx_popup();
    update_selected_stim_chan_popup();
    update_set_textbox();
    
    % Define three callbacks, one for each popup.
    function selected_stimfile_popup_callback()
        % Update the selected_stim_idx_popup string to reflect new choices
        update_selected_stim_idx_popup();
        update_selected_stim_chan_popup();
        update_set_textbox();
        % Call the next popup callback to trigger a redraw
        selected_stim_idx_popup_callback();
    end
    
    function selected_stim_idx_popup_callback()
        global NARFGUI;
        % Call the plot function again via the plot_popup        
        hgfeval(get(NARFGUI{mod_idx}.plot_popup,'Callback'), mod_idx, []);
        drawnow;
    end
    
    function selected_stim_chan_popup_callback()
        global NARFGUI;
        % Call the plot function again via the plot_popup 
        hgfeval(get(NARFGUI{mod_idx}.plot_popup,'Callback'), mod_idx, []);
        drawnow;
    end
end

% This module can be run if all necessary fields have been defined in the
% topmost part of the stack
function isready = isready_baphy(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    isready = (length(stack) == 1) && ...
              all(isfield(x, {'cellid', 'training_set', 'test_set'}));
end

end

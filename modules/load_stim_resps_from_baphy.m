function m = load_stim_resps_from_baphy(args)
% LOAD_STIM_RESPS_FROM_BAPHY
% A module to load stimulus and response files from BAPHY. 
% Returns a function module which implements the MODULE interface.
% See documentation for more information about what a MODULE is.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @load_stim_resps_from_baphy;
m.name = 'load_stim_resps_from_baphy';
m.fn = @do_load_from_baphy;
m.pretty_name = 'Load stim+resp from BAPHY';
m.editable_fields = {'raw_stim_fs', 'raw_resp_fs', 'include_prestim', ...
                     'stimulus_format','stimulus_channel_count', ...
                     'output_stim', 'output_stim_time', ...
                     'output_resp', 'output_resp_time', 'output_respavg'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.raw_stim_fs = 100000;
m.raw_resp_fs = 200;
m.include_prestim = 1;
m.exclude_target_phase = 1;
m.stimulus_channel_count=0; % 0 should be 'autodetect'
m.stimulus_format = 'wav';  % Can be 'wav' or 'envelope'
m.output_stim = 'stim';
m.output_stim_time = 'stim_time';
m.output_resp = 'resp';
m.output_resp_time = 'resp_time';
m.output_respavg = 'respavg';
m.include_prestim = 1;
m.is_data_loader = true; % Special marker used by jackknifing routine

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

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

% ------------------------------------------------------------------------
% Define the 'methods' of this module, as if it were a class

function x = do_load_from_baphy(mdl, x, stack, xxx)
    
    % Merge the training and test set names, which may overlap
    files_to_load = unique({x.training_set{:}, x.test_set{:}});
    
    % Returns a list of raw files with this cell ID
    [cfd, cellids, cellfileids] = dbgetscellfile('cellid', x.cellid);
    x.cfd = cfd;
    
    % If there is not exactly one cell file returned, throw an error.
    if ~isequal(length(cellids), 1)
        error('BAPHY gave %d cellids yet I want only 1.', length(cellids));
    end
    
    % Create the 'dat' cell array and its entries, if it doesn't exist
    len = length(files_to_load);
    %x.dat = cell(1, len);
    if ~isfield(x, 'dat')
        x.dat=struct();
    end
    for f_idx = 1:len;
        f = files_to_load{f_idx};
        
        % SVD 2013-03-08 - check for special est/val subsets
        if strcmpi(f((end-3):end),'_val'),
            datasubset=2;
            fname=f(1:(end-4));
        elseif strcmpi(f((end-3):end),'_est'),
            datasubset=1;
            fname=f(1:(end-4));
        else
            datasubset=0;
            fname=f;
        end
        
        % Find the index number of cfd corresponding to the file
        idx = 0;
        for j = 1:length(cfd)
            if isequal(cfd(j).stimfile, fname)
                idx = j;
            end
        end
        if idx == 0 
            error('Not found in cfd: %s\n', fname); 
        end
        
        % special case kludge to deal with messed up TSP files
        loadbytrial=0;
        a=regexp(fname,'por(\d*).*TSP.*$','tokens');
        if ~isempty(a),
            a=str2num(a{1}{1});
            if a>=31 && a<=65,
                loadbytrial=1;
            end
        end
        
        % Load the raw_stim part of the data structure
        stimfile = [cfd(idx).stimpath cfd(idx).stimfile];
        fprintf('Loading stimulus: %s\n', stimfile);
        if loadbytrial,
            options=[];
            options.filtfmt=mdl.stimulus_format;
            options.fsout=mdl.raw_stim_fs;
            options.chancount=mdl.stimulus_channel_count;
            [stim,stimparam] = loadstimbytrial(stimfile,options);
        else
            [stim,stimparam] = loadstimfrombaphy(stimfile, [], [], ...
                  mdl.stimulus_format, mdl.raw_stim_fs, ...
                  mdl.stimulus_channel_count, 0, mdl.include_prestim);
        end
        
        if strcmp(mdl.stimulus_format, 'wav')
            stim = squeeze(stim);
        end
        
        % TODO: Right now the envelope returned by baphy is not necessarily
        % positive semidefinite, and when run through a square root
        % compressor results in complex numbers being developed. This
        % should be fixed on the baphy side, but for now let's just
        % add a workaround here.
        if strcmp(mdl.stimulus_format, 'envelope') ||...
                strcmpi(mdl.stimulus_format, 'gamma') ||...
                strcmpi(mdl.stimulus_format, 'specgram') ||...
                strcmpi(mdl.stimulus_format, 'specgramv'),
            stim = permute(stim, [2 3 1]);
            stim = abs(stim);
        end
        
        % If you ever need arbitrary dimensions, I started planning
        % out how that could be accomplished...but it doesn't help plotting
        % See documentation in narf/doc/
        % stim = define_dims(stim, {'time', 'stim', 'chan'});
        
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
        
        % Load the raw_resp part of the data structure
        options = [];
        options.includeprestim = mdl.include_prestim;
        options.unit = cfd(idx).unit;
        options.channel  = cfd(idx).channum;
        options.rasterfs = mdl.raw_resp_fs;
        respfile = [cfd(idx).path, cfd(idx).respfile];
        fprintf('Loading response: %s\n', respfile);
        if loadbytrial,
            options.tag_masks={'SPECIAL-TRIAL'};
            [resp, tags,trialset] = loadspikeraster(respfile, options);
            resp=permute(resp,[1 3 2]);
            
            % only keep correct trials
            stim=stim(:,trialset,:);
            
            % mask out responses that aren't in valid stimulus
            % range (ie, responses during & after targets)
            for trialidx=1:size(stim,2),
                gg=find(~isnan(stim(:,trialidx,1)));
                resp((gg(end)+1):end,1,trialidx)=nan;
            end
            
            % find identical trials and call them repeats
            CC=zeros(size(stim,2));
            for dd=1:size(stim,3),
                ts=stim(:,:,dd);
                ts(isnan(ts))=0;
                CC=CC+corrcoef(ts)./size(stim,3);
            end
            CC=CC-triu(CC);
            [t2,t1]=find(CC==1);
            keepidx=1:size(stim,2);
            for ii=unique(t1)',
                mm=t2(find(t1==ii & t2>0));
                tr=squeeze(resp(:,1,mm));
                repcount=size(tr,2)+1;
                if repcount>size(resp,2),
                    resp=cat(2,resp,nan*ones(size(resp,1),...
                             repcount-size(resp,2),size(resp,3)));
                end
                resp(:,2:repcount,ii)=tr;
                keepidx=setdiff(keepidx,mm);
                t2(ismember(t2,mm))=0;
            end
            resp=resp(:,:,keepidx);
            stim=stim(:,keepidx,:);
            
        else
            options.tag_masks={'Reference'};
            [resp, tags] = loadspikeraster(respfile, options);
        end
        
        %keyboard
        
        % SVD pad response with nan's in case reference responses
        % were truncated because of target overlap during behavior.
        % this is a kludge that may be fixed some day in loadspikeraster
        % 2013-05-31: Ivar wraps an if statement around this because this
        % destroyed our ability to sample stim and resp at different rates
        if (mdl.raw_resp_fs == mdl.raw_stim_fs)
            if size(resp,1)<size(stim,1),
                resp((end+1):size(stim,1),:,:)=nan;
            elseif size(resp,1)>size(stim,1),
                resp=resp(1:size(stim,1),:,:);
            end
        end
        
        % SVD 2013-03-08 - if specified, pull out either estimation (fit) or
        % validation (test) subset of the data
        if datasubset,
            repcount=squeeze(sum(~isnan(resp(1,:,:)),2));
            %if max(repcount)>2,
            %    validx=min(find(repcount==max(repcount)));
            %else
                [ff,ii]=sort(repcount,'descend');
                ff=cumsum(ff)./sum(ff);
                validx=ii(1:min(find(ff>=1/15)));
            %end
            %keyboard
            if datasubset==2
                keepidx=validx;
            else
                keepidx=setdiff(find(repcount>0),validx);
            end
            stim=stim(:,keepidx,:);
            resp=resp(:,:,keepidx);
        end
        
        nanstim=find(isnan(stim(:)));
        if ~isempty(nanstim),
            disp('WARNING: nan-valued stims are being set to zero');
            stim(nanstim)=0;
        end
        x.dat.(f).(mdl.output_stim) = stim;
        x.dat.(f).(mdl.output_resp) = permute(resp, [1, 3, 2]);
        x.dat.(f).(mdl.output_respavg) = squeeze(nanmean(x.dat.(f).(mdl.output_resp), ...
                                            3));
                                        
        % Scale respavg so it is a spike rate in Hz
        x.dat.(f).(mdl.output_respavg) = (mdl.raw_resp_fs) .* x.dat.(f).(mdl.output_respavg);
        
        % Create time signals for later convenience
        [s1 s2 s3] = size(x.dat.(f).(mdl.output_stim));
        [r1 r2 r3] = size(x.dat.(f).(mdl.output_resp));
        [a1 a2]    = size(x.dat.(f).(mdl.output_respavg));
        
        % TODO: Check stim, resp, raw_time signal sizes match.       
        x.dat.(f).(mdl.output_stim_time) = (1/mdl.raw_stim_fs).*[1:s1]';
        x.dat.(f).(mdl.output_resp_time) = (1/mdl.raw_resp_fs).*[1:r1]';
    end
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
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);        
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.output_stim_time, mdls{1}.output_stim, ...
            sel, 'Time [s]', what_is_ylabel(mdls));
end

function do_plot_single_stim_channel(sel, stack, xxx)   
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);
    do_plot(xouts, mdls{1}.output_stim_time, mdls{1}.output_stim, ...
            sel, 'Time [s]', what_is_ylabel(mdls)); 
end

function do_plot_respavg(sel, stack, xxx)
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    xouts = xxx(end);
    sel.chan_idx = 1;    
    do_plot(xouts, mdls{1}.output_resp_time, mdls{1}.output_respavg, ...
            sel, 'Time [s]', 'Spike Rate Average [Hz]'); 
end

function do_plot_response_raster(sel, stack, xxx)    
    %[mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
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
        if isfield(x, 'filecodes') && ~isempty(x.filecodes) && ~isempty(x.test_set)
            fc = x.filecodes(or(strcmp(sf, x.training_set), ...
                                strcmp(sf, x.test_set)));
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

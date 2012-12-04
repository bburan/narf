function m = load_stim_resps_from_baphy(args)
% A module to load stimulus and response files from BAPHY. 
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @load_stim_resps_from_baphy;
m.name = 'load_stim_resps_from_baphy';
m.fn = @do_load_from_baphy;
m.pretty_name = 'Load stim+resp from BAPHY';
m.editable_fields = {'raw_stim_fs', 'raw_resp_fs', 'include_prestim'};
m.isready_pred = @module_isready;

% Module fields that are specific to THIS MODULE
m.raw_stim_fs = 100000;
m.raw_resp_fs = 200;
m.include_prestim = 1;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_stim;
m.plot_fns{1}.pretty_name = 'Stimulus vs Time';
m.plot_fns{2}.fn = @do_plot_stim_log_spectrogram;
m.plot_fns{2}.pretty_name = 'Stimulus Log Spectrogram';
m.plot_fns{3}.fn = @do_plot_respavg;
m.plot_fns{3}.pretty_name = 'Response Average vs Time';
m.plot_fns{4}.fn = @do_plot_response_rastered;
m.plot_fns{4}.pretty_name = 'Response Raster Plot';
m.plot_fns{5}.fn = @do_plot_spectro_and_raster;
m.plot_fns{5}.pretty_name = 'Spectrogram + Raster';
m.plot_gui_create_fn = @create_gui;

% ------------------------------------------------------------------------
% Define the 'methods' of this module, as if it were a class
function x = do_load_from_baphy(stack, xxx)
    n = length(stack);
    mdl = stack{n};
    x = xxx{n};
    
    % Merge the training and test set names, which may overlap
    files_to_load = unique({x.training_set{:}, x.test_set{:}});
    
    % Returns a list of raw files with this cell ID
    [cfd, cellids, cellfileids] = dbgetscellfile('cellid', x.cellid);
    x.cfd = cfd;
    
    % If there is not exactly one cell file returned, throw an error.
    if ~isequal(length(cellids), 1)
        error('BAPHY gave %d cellids yet I want only 1.', length(cellids));
    end
    
    % Create the 'dat' cell array and its entries
    len = length(files_to_load);
    x.dat = cell(1, len);
    for f_idx = 1:len;
        f = files_to_load{f_idx};
        
        % Find the index number of cfd corresponding to the file
        idx = 0;
        for j = 1:length(cfd)
            if isequal(cfd(j).stimfile, f)
                idx = j;
            end
        end
        if idx == 0 
            error('Not found in cfd: %s\n', f); 
        end        
          
        % Load the raw_stim part of the data structure
        stimfile = [cfd(idx).stimpath cfd(idx).stimfile];
        fprintf('Loading stimulus: %s\n', stimfile);
        stim = loadstimfrombaphy(stimfile, [], [], 'wav', ...
            mdl.raw_stim_fs, 1, 0, mdl.include_prestim);
        [d1 d2 d3] = size(stim);
        if d1 ~= 1
            log_dbg('Stimulus matrix was: [%d %d %d]\n', d1, d2, d3);
            log_err('Stimulus size was not [1xNxS]!?\n');
        end
        x.dat.(f).raw_stim = permute(stim, [3 2 1]);
        x.dat.(f).raw_stim_fs = mdl.raw_stim_fs;         % TODO: Remove this SPOT violation
        x.dat.(f).include_prestim = mdl.include_prestim; % TODO: Remove this SPOT violation
        
        % Load the raw_resp part of the data structure
        options = [];
        options.includeprestim = mdl.include_prestim;
        options.unit = cfd(idx).unit;
        options.channel  = cfd(idx).channum;
        options.rasterfs = mdl.raw_resp_fs;
        respfile = [cfd(idx).path, cfd(idx).respfile];
        fprintf('Loading response: %s\n', respfile);
        [resp, tags] = loadspikeraster(respfile, options);
        x.dat.(f).raw_resp = permute(resp, [3 1 2]);
        x.dat.(f).raw_respavg = squeeze(sum(x.dat.(f).raw_resp, 3));
        x.dat.(f).raw_resp_fs = mdl.raw_resp_fs;
        
        % Check raw_stim, raw_resp, and raw_time signal sizes match.
        % TODO: This is no longer relevant if stim and resp are sampled at
        % different frequencies!
%         if ~(isequal(s1, r1, a1) & isequal(s2, r2, a2))
%             fprintf('Stim [%d %d], Resp: [%d %d %d] Respavg=[%d %d]\n',...
%                 s1,s2,r1,r2,r3, a1,a2);
%             error('Stimulus, Response, and Average matrices size mismatch.\n');
%         end
        
        % Create time signals for convenience
        [s1 s2]    = size(x.dat.(f).raw_stim);
        [r1 r2 r3] = size(x.dat.(f).raw_resp);
        [a1 a2]    = size(x.dat.(f).raw_respavg);
        x.dat.(f).raw_stim_time = (1/mdl.raw_stim_fs).*[1:s2]';
        x.dat.(f).raw_resp_time = (1/mdl.raw_resp_fs).*[1:r2]';
        
    end
end

function do_plot_stim(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Read the GUI to find out the selected stim files
    c = cellstr(get(mdl.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(mdl.plot_gui.selected_stimfile_popup, 'Value')};
    idx = get(mdl.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);

    % TODO: If there was any problem in the above, report it
    
    cla;
    plot(dat.raw_stim_time, dat.raw_stim(idx,:), 'k-');
    axis tight;    
end

function do_plot_stim_log_spectrogram(stack, xxx)
    mdl = stack{end};    
    x = xxx{end};
    
    % Read the GUI to find out the selected stim files
    c = cellstr(get(mdl.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(mdl.plot_gui.selected_stimfile_popup, 'Value')};
    idx = get(mdl.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
    logfsgram(dat.raw_stim(idx,:)', 4048, mdl.raw_stim_fs, [], [], 500, 12); 
    caxis([-20,40]);  % TODO: use a 'smarter' caxis here
    axis tight;
    
end

function do_plot_respavg(stack, xxx)
    mdl = stack{end};    
    x = xxx{end};
    
    % Read the GUI to find out the selected stim files
    c = cellstr(get(mdl.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(mdl.plot_gui.selected_stimfile_popup, 'Value')};
    idx = get(mdl.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
	
    plot(dat.raw_resp_time(:), dat.raw_respavg(idx,:), 'k-');
    %[xs,ys] = find(dat.raw_respavg(idx, :) > 0);
    % bar(dat.raw_resp_time(ys), dat.raw_respavg(idx,ys), 0.01, 'k-');
    % axis([0 dat.raw_resp_time(end) 0 2]);
    axis([0 dat.raw_resp_time(end) 0 max(dat.raw_respavg(idx,:))]);
end

function do_plot_response_rastered(stack, xxx)
    mdl = stack{end};    
    x = xxx{end};
    
    % Read the GUI to find out the selected stim files
    c = cellstr(get(mdl.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(mdl.plot_gui.selected_stimfile_popup, 'Value')};
    idx = get(mdl.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    [S, N, R] = size(dat.raw_resp);
    cla;
    hold on;
    for j = 1:R
        [xs,ys] = find(dat.raw_resp(idx, :, j) > 0);
        plot(dat.raw_resp_time(ys), j*dat.raw_resp(idx,ys,j), 'k.');
	end
    axis([0 dat.raw_resp_time(end) 0 R+1]);
    %setAxisLabelCallback(gca, @(y)(y), 'Y');
    hold off;
end

function do_plot_spectro_and_raster(stack, xxx)
    mdl = stack{end};    
    x = xxx{end};
    
    % Read the GUI to find out the selected stim files
    c = cellstr(get(mdl.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(mdl.plot_gui.selected_stimfile_popup, 'Value')};
    idx = get(mdl.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    cla;
    hold on;
    % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
    logfsgram(dat.raw_stim(idx,:)', 4048, mdl.raw_stim_fs, [], [], 500, 12); 
    caxis([-20,40]);  % TODO: use a 'smarter' caxis here
    h = get(gca, 'YLim');
    d = h(2) - h(1);
    axis tight;
    [S, N, R] = size(dat.raw_resp);
    for j = 1:R
        [xs,ys] = find(dat.raw_resp(idx, :, j) > 0);
        plot(dat.raw_resp_time(ys), ...
             h(1) + (j/R)*d*(d/(d+d/R))*dat.raw_resp(idx,ys,j), 'k.');
    end
    hold off;
end


function hs = create_gui(parent_handle, stack, xxx)
    pos = get(parent_handle, 'Position');
    w = pos(3) - 10;
    h = pos(4) - 10;
    hs = [];

    %mod_idx = length(stack);
    m = stack{end};
    mod_idx = length(stack);
    x = xxx{end};
    stimfiles = fieldnames(x.dat);
    
    % Create a popup which selects
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'StimFile:', ...
        'Units', 'pixels', 'Position', [5 (h-25) w 25]);
    hs.selected_stimfile_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [5 (h-50) w 25], ...
        'Callback', @(a,b,c) selected_stimfile_popup_callback());
    
    % Create a stimfile selector
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left', 'String', 'Stim Num:', ...
        'Units', 'pixels', 'Position', [5 (h-75) w 25]);
    hs.selected_stim_idx_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'NONE', ...
        'Units', 'pixels', 'Position', [5 (h-100) w 25], ...
        'Callback', @(a,b,c) selected_stim_idx_popup_callback());

    % Two functions to populate the two popup menus
    function update_selected_stimfile_popup()
        fns = fieldnames(x.dat);
        set(hs.selected_stimfile_popup, 'String', char(fns));
    end
    
    function update_selected_stim_idx_popup()
        % Get the selected stim files
        c = cellstr(get(hs.selected_stimfile_popup, 'String'));
        sf = c{get(hs.selected_stimfile_popup, 'Value')};
        
        if isfield(x.dat, sf)
            [d1, d2] = size(x.dat.(sf).raw_respavg);
            d = {};
            for i = 1:d1
                d{i} = sprintf('%d',i);
            end
            set(hs.selected_stim_idx_popup, 'String', char(d));
            set(hs.selected_stim_idx_popup, 'Value', 1);
        else
            error('Selected stimulus file not found: %s', sf);
        end
    end
    
    % Call the two update functions once to build the lists
    update_selected_stimfile_popup();
    update_selected_stim_idx_popup();
    
    % Define two callbacks, one for each popup.
    function selected_stimfile_popup_callback()
        % Update the selected_stim_idx_popup string to reflect new choices
        update_selected_stim_idx_popup();
        
        % Call the next popup callback to trigger a redraw
        selected_stim_idx_popup_callback();
    end
    
    function selected_stim_idx_popup_callback()
        % Call the plot function again via a sneaky, undocumented callback
        hgfeval(get(m.gh.plot_popup,'Callback'), mod_idx, []);
        drawnow;
    end
end

% This module can be run if all necessary fields have been defined in the
% topmost part of the stack
function isready = module_isready(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    isready = (length(stack) == 1) && ...
              all(isfield(x, {'cellid', 'training_set', 'test_set'}));
end

end
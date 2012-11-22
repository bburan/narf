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
m.plot_fns = {'Special Plot', @do_plot_special};
m.isready_pred = @module_isready;

% Optional fields
m.plot_gui_create_fn = @create_gui;
m.plot_gui = [];

% Module fields that are specific to THIS MODULE
m.raw_stim_fs = 100000;
m.raw_resp_fs = 10000;
m.include_prestim = 1;
m.training_set = {};
m.t_set = {};

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Values computed below this point are not directly editable. Instead, 
% they are computed from the above values in a deterministic manner.

% NONE for this module.

% Finally, define the 'methods' of this module, as if it were a class
function x = do_load_from_baphy(stack, x)
    mdl = stack{end};
    
    % Merge the training and test set names, which may overlap
    files_to_load = unique({mdl.training_set{:}, mdl.test_set{:}});
    
    % Create the 'dat' cell array and its entries
    len = length(files_to_load);
    x.dat = cell(1, len);
    for f_idx = 1:len;
        f = files_to_load{f_idx};
        
        % Find the index number of mdl.cfd corresponding to the file
        idx = 0;
        for j = 1:length(mdl.cfd)
            if isequal(mdl.cfd(j).stimfile, f)
                idx = j;
            end
        end
        if idx == 0 
            error('Not found in mdl.cfd: %s', f); 
        end        
          
        % Load the raw_stim part of the data structure
        stimfile = [mdl.cfd(idx).stimpath mdl.cfd(idx).stimfile];
        log_msg('Loading stimulus: %s', stimfile);
        stim = loadstimfrombaphy(stimfile, [], [], 'wav', ...
            mdl.raw_stim_fs, 1, 0, mdl.include_prestim);
        [d1 d2 d3] = size(stim);
        if d1 ~= 1
            log_dbg('Stimulus matrix was: [%d %d %d]', d1, d2, d3);
            log_err('Stimulus size was not [1xNxS]!?');
        end
        x.dat.(f).raw_stim = permute(stim, [3 2 1]);
        x.dat.(f).raw_stim_fs = mdl.raw_stim_fs;
        x.dat.(f).include_prestim = mdl.include_prestim;
        
        % Load the raw_resp part of the data structure
        options = [];
        options.includeprestim = includeprestim;
        options.unit = mdl.cfd(idx).unit;
        options.channel  = mdl.cfd(idx).channum;
        options.rasterfs = mdl.raw_resp_fs;
        respfile = [mdl.cfd(idx).path, mdl.cfd(idx).respfile];
        log_msg('Loading response: %s', respfile);
        [resp, tags] = loadspikeraster(respfile, options);
        x.dat.(f).raw_resp = permute(resp, [3 1 2]);
        x.dat.(f).raw_respavg = squeeze(sum(x.dat.(f).raw_resp, 3));
        x.dat.(f).raw_resp_fs = mdl.raw_resp_fs;
        
        % Check raw_stim, raw_resp, and raw_time signal sizes match.
        [s1 s2]    = size(x.dat.(f).raw_stim);
        [r1 r2 r3] = size(x.dat.(f).raw_resp);
        [a1 a2]    = size(x.dat.(f).raw_respavg);
        if ~(isequal(s1, r1, a1) & isequal(s2, r2, a2))
            log_dbg('Stim [%d %d], Resp: [%d %d %d] Respavg=[%d %d]',...
                s1,s2,r1,r2,r3, a1,a2);
            log_err('Stimulus, Response, and Average matrices size mismatch.');
        end
        
        % Create time signals for convenience
        x.dat.(f).raw_stim_time = (1/mdl.raw_stim_fs).*[1:s2]';
        x.dat.(f).raw_resp_time = (1/mdl.raw_resp_fs).*[1:r2]';
        
    end
end

function do_plot_stim(stack, x)
    mdl = stack{end};
    dat = x.dat.(mdl.plot_gui.selected_stimfile);
    cla;
    plot(dat.raw_stim_time, dat.raw_stim(mdl.plot_gui.selected_stim_idx,:), 'k-');
    axis tight;    
end

function do_plot_stim_log_spectrogram(stack, x)
    mdl = stack{end};    
    dat = x.dat.(mdl.plot_gui.selected_stimfile);
    cla;
    % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
    logfsgram(dat.raw_stim(mdl.plot_gui.stim_idx,:)', 4048, 100000, [], [], 500, 12); 
    % TODO: Remove 4048, 100000 here and use global
    caxis([-20,40]);  % TODO: use a 'smarter' caxis here
    axis tight;
    
end

function do_plot_respavg(stack, x)
    mdl = stack{end};
    dat = x.dat.(mdl.plot_gui.selected_stimfile);
    stim_idx = mdl.plot_gui.selected_stim_idx;
    [S, N, R] = size(dat.raw_resp);
	[xs,ys] = find(dat.raw_respavg(stim_idx, :) > 0);
    bar(dat.raw_time(ys), dat.raw_respavg(stim_idx,ys), 0.01, 'k');
    axis([0 dat.raw_time(end) 0 2]);
end

function do_plot_response_rastered(stack, x)
    mdl = stack{end};
    dat = x.dat.(mdl.plot_gui.selected_stimfile);
    stim_idx = mdl.plot_gui.selected_stim_idx;
    [S, N, R] = size(dat.raw_resp);
    axes(handles.resp_view_axes); cla; hold on;
    for j = 1:R
        [xs,ys] = find(dat.raw_resp(stim_idx, :, j) > 0);
        plot(dat.raw_time(ys), j/R*dat.raw_resp(stim_idx,ys,j), 'k.');
	end
    axis([0 dat.raw_time(end) 1/(2*R) 1+(1/(2*R))]);
end

function hs = create_gui(parent_handle)
    pos = get(parent_handle, 'Position');
    w = pos(3);
    h = pos(4);

    % Create a popup which selects
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'String', 'StimFile:', 'Units', 'pixels', 'Position', [5 (h-20) 50 25]);
    select_stimfile_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'NONE', ...
        'Units', 'pixels', 'Position', [55 (h-20) 80 25], ...
        'Callback', @(a,b,c) disp('FN popup callback TODO.'));
    
    % Create a stimfile selector
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'String', 'Stim Num:', 'Units', 'pixels', 'Position', [5 (h-60) 50 25]);
    select_stim_idx_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'NONE', ...
        'Units', 'pixels', 'Position', [55 (h-60) 80 25], ...
        'Callback', @(a,b,c) disp('FN popup callback TODO.'));
    
    hs = [];
    hs.selected_stimfile = select_stimfile_popup;
    hs.selected_stim_idx = select_stim_idx_popup;

end

% This module can be run if all necessary fields have been defined in the
% topmost part of the stack
function isready = module_isready(stack, x)
    mdl = stack{end};
    isready = all(isfield(mdl, {'cellid', 'training_set', 'test_set'}));
end

end
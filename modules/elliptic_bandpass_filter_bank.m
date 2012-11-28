function m = elliptic_bandpass_filter_bank(args)
% An elliptic bandpass filter bank module creation function.
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @elliptic_bandpass_filter_bank;
m.name = 'elliptic_bandpass_filter_bank';
m.fn = @do_elliptic_filter;
m.pretty_name = 'Elliptic Bandpass Filter Bank';
m.editable_fields = {'low_freqs', 'high_freqs', 'order', 'sampfs', 'stop_dB'};
m.isready_pred = @preproc_filter_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim;
m.plot_fns{1}.pretty_name = 'Filtered Stimulus vs Time';
m.plot_fns{2}.fn = @do_plot_filtered_spectrogram;
m.plot_fns{2}.pretty_name = 'Filtered Stimulus Spectrogram';
m.plot_fns{3}.fn = @do_plot_elliptic_bandpass_filter_bank_frq_resp;
m.plot_fns{3}.pretty_name = 'Filter Frequency Response';
m.plot_gui_create_fn = @create_gui;

% Module fields that are specific to THIS MODULE
m.low_freqs = [2000 20000];  % Bottom frequencies of bandpass filters
m.high_freqs = [4000 27000]; % Top frequencies of bandpass filters
m.order = 4;                 % What order should the filter be?
m.sampfs = 100000;           % TODO: REMOVE THIS AND AUTODETECT IT
m.stop_dB = 50;              % Ratio of passband/stopband attenuation

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Values computed from here on are not directly editable, but are based on
% the above values in a deterministic manner.
m.coefs={}; 
for i = 1:length(m.low_freqs);
    [B,A] = ellip(m.order, 0.5, m.stop_dB, ...
                  [(m.low_freqs(i)/m.sampfs)*2,...
                   (m.high_freqs(i)/m.sampfs)*2]);   
    m.coefs{i} = {B,A};
end 

% Finally, define the 'methods' of this module, as if it were a class
function x = do_elliptic_filter(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % For each data file...
    sfs = fieldnames(x.dat);
    for s_idx = 1:length(sfs)
        % make a matrix to store all the filter responses...
        filtered_x = [];
        for idx = 1:length(mdl.low_freqs)
            tmp = filter(mdl.coefs{idx}{1}, mdl.coefs{idx}{2}, x.dat.(sfs{s_idx}).raw_stim,[],2);
            filtered_x = cat(3, filtered_x, tmp); 
        end
        % Store that matrix in our data structure
        x.dat.(sfs{s_idx}).pp_stim = filtered_x;
    end
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
    dat = x.dat.(sf);
    filt_idx = get(mdl.plot_gui.selected_elliptic_filter_popup, 'Value');
    
    plot(dat.raw_stim_time, ...
         squeeze(dat.pp_stim(stim_idx,:,filt_idx)), ...
         pickcolor(filt_idx));
    axis tight;
    drawnow;
end

function do_plot_filtered_spectrogram(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    filt_idx = get(mdl.plot_gui.selected_elliptic_filter_popup, 'Value');
    
    logfsgram(dat.pp_stim(stim_idx,:, filt_idx)', 4048, baphy_mod.raw_stim_fs, [], [], 500, 12); 
    caxis([-20,40]);
    drawnow;
end
        
function do_plot_elliptic_bandpass_filter_bank_frq_resp(stack, xxx)   
    mdl = stack{end};
    x = xxx{end};

    hold on;
    for filt_idx = 1:length(mdl.low_freqs)
        ww = 0:(pi/1000):pi;
        H = freqz(mdl.coefs{filt_idx}{1}, mdl.coefs{filt_idx}{2}, ww);
        loglog(ww, abs(H), pickcolor(filt_idx));
        setAxisLabelCallback(gca, @(f) (f*mdl.sampfs/(3.14*2)), 'X');
        axis tight;
    end 
    hold off;
end


function hs = create_gui(parent_handle, stack, xxx)
    pos = get(parent_handle, 'Position');
    w = pos(3) - 10;
    h = pos(4) - 10;
    hs = [];

    mdl = stack{end};
    mod_idx = length(stack);
    x = xxx{end};
    
    % Create a popup which selects
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'Filter#:', ...
        'Units', 'pixels', 'Position', [5 (h-25) w 25]);
    hs.selected_elliptic_filter_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [5 (h-50) w 25], ...
        'Callback', @(a,b,c) selected_elliptic_filter_popup_callback());
        
    % Fill that popup with the number of filters
    d = {};
    for ii = 1:length(mdl.low_freqs)
        d{ii} = sprintf('%d',ii);
    end
    set(hs.selected_elliptic_filter_popup, 'String', char(d));
    set(hs.selected_elliptic_filter_popup, 'Value', 1);
    
    function selected_elliptic_filter_popup_callback()
        % Call the plot function again via a sneaky, undocumented callback
        hgfeval(get(mdl.gh.plot_popup,'Callback'), mod_idx, []);
        drawnow;
    end
end

end
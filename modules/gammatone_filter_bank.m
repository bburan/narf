function m = gammatone_filter_bank(args)
% A gammatone filter bank module creation function.
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @gammatone_filter_bank;
m.name = 'gammatone_filter_bank';
m.fn = @do_gammatone_filter;
m.pretty_name = 'Gammatone Filter Bank';
m.editable_fields = {'bank_min_freq', 'bank_max_freq', 'num_gammatone_filters', 'align_phase'};
m.isready_pred = @preproc_filter_isready;

% Module fields that are specific to THIS MODULE
m.bank_min_freq = 100;
m.bank_max_freq = 30000;
m.num_gammatone_filters = 32;
m.align_phase = 0;
m.raw_stim_freq = 100000;
           
% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim;
m.plot_fns{1}.pretty_name = 'Filtered Stimulus vs Time';
m.plot_fns{2}.fn = @do_plot_filtered_spectrogram;
m.plot_fns{2}.pretty_name = 'Filtered Stimulus Spectrogram';
m.plot_fns{3}.fn = @do_plot_frequency_response;
m.plot_fns{3}.pretty_name = 'Filter Frequency Responses';
m.plot_fns{4}.fn = @do_plot_frequency_response;
m.plot_fns{4}.pretty_name = 'Gammatonegram';
m.plot_gui_create_fn = @create_gui;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Finally, define the 'methods' of this module, as if it were a class
function x = do_gammatone_filter(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    sfs = fieldnames(x.dat);
    for s_idx = 1:length(sfs)
        [gamma_resp, gamma_envs, gamma_frqs] = ...
            gammatonebank(x.dat.(sfs{s_idx}).raw_stim(s_idx, :), ...
                          m.bank_min_freq, m.bank_max_freq, ...
                          m.num_gammatone_filters, baphy_mod.raw_stim_fs, ...
                          m.align_phase);                  
        x.dat.(sfs{s_idx}).pp_stim = gamma_resp;
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

function do_plot_frequency_response(stack, xxx)   
    mdl = stack{end};
    x = xxx{end};
    % TODO
    plot([1,2,3],[3,4,3]);
end

function do_plot_gammatone_filter_as_colormap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find log intensity and smooth slightly
    %     gamma_sqr = [gamma_resp.^2];        
    %     gamma_pow = zeros(N_gfs, N_cols);
    %     for i = 1:N_cols
    %         gamma_pow(:,i) = 20*log10(sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2)));
    %         %gamma_pow(:,i) = sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2));
    %     end
    % TODO:
    % imagesc(gamma_pow);
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
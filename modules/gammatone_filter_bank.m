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
m.editable_fields = {'bank_min_freq', 'bank_max_freq', 'num_channels', 'align_phase'};
m.isready_pred = @preproc_filter_isready;

% Module fields that are specific to THIS MODULE
m.bank_min_freq = 500;
m.bank_max_freq = 30000;
m.num_channels = 5;
m.align_phase = false;
m.raw_stim_freq = 100000;
           
% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim; % TODO: make this call a generic like plot_time_series('raw_stim_time', 'pp_stim')
m.plot_fns{1}.pretty_name = 'Filtered Stimulus vs Time';
m.plot_fns{2}.fn = @do_plot_filtered_spectrogram;
m.plot_fns{2}.pretty_name = 'Filtered Stimulus Spectrogram';
m.plot_fns{3}.fn = @do_plot_frequency_response;
m.plot_fns{3}.pretty_name = 'Filter Frequency Responses';
m.plot_fns{4}.fn = @do_plot_gammatone_filter_as_colormap;
m.plot_fns{4}.pretty_name = 'Gammatonegram';
m.plot_gui_create_fn = @(h, stk, xx) create_filter_selector_gui(h, stk, xx, m.num_channels);

% Finally, define the 'methods' of this module, as if it were a class
function x = do_gammatone_filter(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    % Exotic way to loop over field names using ' and {1}...
    for sf = fieldnames(x.dat)', sf = sf{1};
        [S, N] = size(x.dat.(sf).raw_stim);
        ret = zeros(m.num_channels, N, S);
        for s = 1:S
            fprintf('%d\n', s);
            [ret(:,:, s), gamma_envs, gamma_frqs] = ...
                gammatonebank(x.dat.(sf).raw_stim(s, :), ...
                              m.bank_min_freq, m.bank_max_freq, ...
                              m.num_channels, baphy_mod.raw_stim_fs, ...
                              m.align_phase);         
            % ret = cat(3, ret, gamma_resp);
        end
        x.dat.(sf).pp_stim = permute(ret, [3,2,1]); 
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
    filt_idx = get(mdl.plot_gui.selected_filter_popup, 'Value');
    
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
    
    filt_idx = get(mdl.plot_gui.selected_filter_popup, 'Value');
    
    logfsgram(dat.pp_stim(stim_idx,:, filt_idx)', 4048, baphy_mod.raw_stim_fs, [], [], 500, 12); 
    caxis([-20,40]);
    drawnow;
end

function do_plot_frequency_response(stack, xxx)   
    mdl = stack{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    % Stupid approximation using white noise
    sr=baphy_mod.raw_stim_fs;   % Sample rate
    noise = wgn(5*sr, 1,0);  % 5 secs of noisy samples

    filt_idx = get(mdl.plot_gui.selected_filter_popup, 'Value');
    frqs = MakeErbCFs(mdl.bank_min_freq, mdl.bank_max_freq, mdl.num_channels);   
    fc = frqs(filt_idx);  % Center frequency of gamma filter
        
    y1 = gammatone(noise, sr, fc, mdl.align_phase);
    L=length(y1);
    NFFT = 2^nextpow2(L);
    Y0 = fft(y1,NFFT);
    Y1 = sqrt((Y0(1:NFFT/2+1)/L).^2);
    faxis = sr/2*linspace(0,1,NFFT/2+1);
    
    % We could plot faxis and Y1 now, but we'll "smooth" first...
    % Just for visualization!
    n_smooth = 50; % How many samples to bring together to smooth. EVEN
    n = length(Y1);
    n_bins = fix(n/n_smooth);
    Y_smoothed = mean(reshape(Y1(1:n_bins*n_smooth), n_smooth, n_bins));
    F_smoothed = mean(reshape(faxis(1:n_bins*n_smooth), n_smooth, n_bins));
    
    % Plot single-sided amplitude spectrum
    loglog(F_smoothed, 20*Y_smoothed, 'k-');
    axis([mdl.bank_min_freq, mdl.bank_max_freq, 10^-6, 1]);
    %loglog(faxis, Y1, 'k-');
    
end

function do_plot_gammatone_filter_as_colormap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    cla;
    % TODO:
    % Find log intensity and smooth slightly
    %     gamma_sqr = [gamma_resp.^2];        
    %     gamma_pow = zeros(N_gfs, N_cols);
    %     for i = 1:N_cols
    %         gamma_pow(:,i) = 20*log10(sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2)));
    %         %gamma_pow(:,i) = sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2));
    %     end
    % imagesc(gamma_pow);
end

end
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
m.editable_fields = {'bank_min_fs', 'bank_max_fs', ...
                     'num_channels', 'align_phase', 'use_env', ...
                     'input', 'time', 'output' };
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.bank_min_fs = 500;
m.bank_max_fs = 30000;
m.num_channels = 10;
m.align_phase = false;
m.use_env = false;
m.raw_stim_fs = 100000;
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim; 
 m.plot_fns{1}.pretty_name = 'Filtered Stimulus vs Time';
% m.plot_fns{2}.fn = @do_plot_gammatone_filter_as_colormap;
% m.plot_fns{2}.pretty_name = 'Gammatonegram';
% m.plot_fns{3}.fn = @do_plot_filtered_spectrogram;
% m.plot_fns{3}.pretty_name = 'Filtered Stimulus Spectrogram';
% m.plot_fns{4}.fn = @do_plot_frequency_response;
% m.plot_fns{4}.pretty_name = 'Filter Frequency Responses';

m.plot_gui_create_fn = @(h, stk, xx) create_filter_selector_gui(h, stk, xx, m.num_channels);

function x = do_gammatone_filter(mdl, x, stack, xxx)

    for sf = fieldnames(x.dat)', sf = sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        ret = zeros(T, C * mdl.num_channels, S);
        for s = 1:S
            %fprintf('gammatone_filter_bank.m: [%d/%d]\n', s, S)
            for c = 1:C
                [gamma_bms, gamma_envs, ~] = ...
                    gammatonebank(x.dat.(sf).(mdl.input)(:, s, c), ...
                                  mdl.bank_min_fs, mdl.bank_max_fs, ...
                                  mdl.num_channels, mdl.raw_stim_fs, ...
                                  mdl.align_phase);    
                b = mdl.num_channels*(c-1) + 1;               
                d = mdl.num_channels*(c);
                if (mdl.use_env)
                    ret(:, b:d, s) = gamma_envs';
                else
                    ret(:, b:d, s) = gamma_bms';
                end
            end
        end
        x.dat.(sf).(mdl.output) = abs(permute(ret, [1,3,2])); 
    end
end

function do_plot_filtered_stim(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));  
    if mdls{1}.use_env
        ylab = 'Basilar Envelope';
    else
        ylab = 'Basilar Motion';
    end    
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', ylab);           
end

% function do_plot_filtered_spectrogram(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     
%     [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);
%     
%     filt_idx = get(mdl.plot_gui.selected_filter_popup, 'Value');
%     
%     logfsgram(dat.(mdl.output)(stim_idx,:, filt_idx)', 4048, baphy_mod.raw_stim_fs, [], [], 500, 12); 
%     caxis([-20,40]);
% end
% 
% function do_plot_frequency_response(stack, xxx)   
%     mdl = stack{end};
%     
%     [baphy_mod, ~] = find_modules(stack, 'load_stim_resps_from_baphy', true);
%     
%     % Stupid approximation using white noise
%     sr=baphy_mod.raw_stim_fs;   % Sample rate
%     noise = wgn(5*sr, 1,0);  % 5 secs of noisy samples
% 
%     filt_idx = get(mdl.plot_gui.selected_filter_popup, 'Value');
%     frqs = MakeErbCFs(mdl.bank_min_fs, mdl.bank_max_fs, mdl.num_channels);   
%     fc = frqs(filt_idx);  % Center frequency of gamma filter
%         
%     y1 = gammatone(noise, sr, fc, mdl.align_phase);
%     L=length(y1);
%     NFFT = 2^nextpow2(L);
%     Y0 = fft(y1,NFFT);
%     Y1 = sqrt((Y0(1:NFFT/2+1)/L).^2);
%     faxis = sr/2*linspace(0,1,NFFT/2+1);
%     
%     % We could plot faxis and Y1 now, but we'll "smooth" first...
%     % Just for visualization!
%     n_smooth = 50; % How many samples to bring together to smooth. EVEN
%     n = length(Y1);
%     n_bins = fix(n/n_smooth);
%     Y_smoothed = mean(reshape(Y1(1:n_bins*n_smooth), n_smooth, n_bins));
%     F_smoothed = mean(reshape(faxis(1:n_bins*n_smooth), n_smooth, n_bins));
%     
%     % Plot single-sided amplitude spectrum
%     loglog(F_smoothed, 20*Y_smoothed, 'k-');
%     axis([mdl.bank_min_fs, mdl.bank_max_fs, 10^-6, 1]);
%     %loglog(faxis, Y1, 'k-');
%     
% end
% 
% function do_plot_gammatone_filter_as_colormap(stack, xxx)
%     
%     [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));  
%              
%     ii=1;  % assuming just a single data set for now...
%     h = imagesc(xouts{ii}.dat.(sel.stimfile).(mdls{1}.time)(:),...
%                 1:size(xouts{ii}.dat.(sel.stimfile).(mdls{1}.output),3),...
%                 squeeze(xouts{ii}.dat.(sel.stimfile).(mdls{1}.output)(:, sel.stim_idx, :))');
%     do_xlabel('Time [s]');
%     do_ylabel('Output [-]');
%     
%     set(gca,'YDir','normal');
%     axis xy tight;
%     % TODO:
%     % Find log intensity and smooth slightly
%     %     gamma_sqr = [gamma_resp.^2];        
%     %     gamma_pow = zeros(N_gfs, N_cols);
%     %     for i = 1:N_cols
%     %         gamma_pow(:,i) = 20*log10(sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2)));
%     %         %gamma_pow(:,i) = sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2));
%     %     end
%     % imagesc(gamma_pow);
% end

end

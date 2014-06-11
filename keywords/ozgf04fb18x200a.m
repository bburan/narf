function ozgf04fb18x200a ()
% One-zero gamma-tone filterbank, 4th order, 18 channels at 200Hz.
% Appends 4 modules to stack:
%    1. load_stim_resps_from_baphy
%    2. pz_wavelet 
%    3. downsample_signal
%    4. normalize_channels

global MODULES STACK XXX META;

load_fs = 50000;
samp_fs = 200;

n_chans = 18;
CFs = logspace(log10(0.2), log10(20), n_chans);
Qs = ones(1, n_chans) * (3*CFs(1) ./ (CFs(2) - CFs(1)));

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', samp_fs, ...
                                  'raw_stim_fs', load_fs,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'wav'))); 
                            
append_module(MODULES.pz_wavelet_filterbank.mdl(...
                struct('zeros', zeros(1,n_chans), ...
                       'N_order', 4, ...
                       'M_order', 1, ...
                       'center_freq_khz', CFs, ...
                       'Q_factor', Qs, ...
                       'delayms', zeros(1,n_chans), ...
                       'align_peak', true, ...
                       'input', 'stim', ...
                       'output', 'stim')));

append_module(MODULES.downsample_signal.mdl(...
        struct('input', 'stim', ...
               'input_freq', load_fs, ...
               'input_time', 'stim_time', ...
               'output', 'stim', ...
               'output_freq', samp_fs, ...
               'output_time', 'stim_time', ...
               'downsampler', @decimate))); % Try @decimate, or @conv_fn, or @resample

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

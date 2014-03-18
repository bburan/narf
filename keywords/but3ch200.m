function but3ch200()

global MODULES STACK;

load_fs = 50000;
samp_fs = 200;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', samp_fs, ...
                                  'raw_stim_fs', load_fs,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'wav'))); 
                              % Module fields that are specific to THIS MODULE

append_module(MODULES.bandpass_filter_bank.mdl(...
        struct('low_freqs', [1], ...
               'high_freqs', [10], ...
               'order', 3, ...
               'sampfs', load_fs, ...
               'stop_dB', 50, ...
               'function', @butter, ...
               'fit_fields', {{'low_freqs', 'high_freqs'}})));
           
append_module(MODULES.downsample_signal.mdl(...
        struct('input', 'stim', ...
               'input_freq', load_fs, ...
               'input_time', 'stim_time', ...
               'output', 'stim', ...
               'output_freq', samp_fs, ...
               'output_time', 'stim_time', ...
               'downsampler', @decimate))); % Try @decimate, or @conv_fn, or @resample

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

fitell();
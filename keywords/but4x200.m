function but4x200()

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
               'order', 4, ...
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

append_module(MODULES.weight_channels.mdl(struct('weights', [1], ...
                                                 'y_offset', 0, ...
                                                 'fit_fields', {{'y_offset', 'weights'}})));

                                             
nmse();
for ii = 1:6  
    fit_scaat('StopAtAbsScoreDelta', 10^-ii, ...          
              'InitStepSize', 10.0, ...
              'StopAtStepsize', 10^-4);
end
pop_module(); % Remove NMSE
pop_module(); % Remove wc01

% Stop fitting the PZ wavelet. 
STACK{end-2}{1}.fit_fields = {};

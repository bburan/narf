function apgt04x200delay ()

global MODULES STACK XXX META;

load_fs = 50000;
samp_fs = 200;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', samp_fs, ...
                                  'raw_stim_fs', load_fs,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'wav'))); 

append_module(MODULES.pz_wavelet.mdl(...
                struct('fit_fields', {{'center_freq_khz', 'Q_factor', 'delayms'}}, ...
                       'N_order', 4, ...
                       'center_freq_khz', 1.7, ...                       
                       'Q_factor', 0.7, ... 
                       'delayms', 4, ... 
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

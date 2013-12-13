function apgt01 ()

global MODULES STACK XXX META;

append_module(MODULES.pz_wavelet.mdl(...
                struct('fit_fields', {{'center_freq_khz', 'Q_factor', ...
                                       'N_order', 'delayms'}}, ...
                       'center_freq_khz', 1.7, ...
                       'N_order', 3, ...
                       'Q_factor', 0.7, ... 
                       'delayms', 4, ... 
                       'input', 'stim', ...
                       'output', 'stim')));  

append_module(MODULES.downsample_signal.mdl(...
                struct('input', 'stim', ...
                       'input_freq', 30000, ...
                       'output', 'stim', ...
                       'output_freq', 200, ...
                       'output_time', 'stim_time')));            

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.weight_channels.mdl(struct('weights', [10.3], ...
                                                 'y_offset', 8.3, ...
                                                 'fit_fields', {{'y_offset', 'weights'}})));

nmse();
fit_boo('StopAtAbsScoreDelta', 10^-3, 'StepGrowth', 1.3);
%keyboard;
[a,b,c,d] = fit_scaat('InitStepSize', 1.0, 'StopAtAbsScoreDelta', 10^-1);
[a,b,c,d] = fit_scaat('InitStepSize', 1.0, 'StopAtAbsScoreDelta', 10^-2);
[a,b,c,d] = fit_scaat('InitStepSize', 1.0, 'StopAtAbsScoreDelta', 10^-3);
[a,b,c,d] = fit_scaat('InitStepSize', 1.0, 'StopAtAbsScoreDelta', 10^-4);

% Remove MSE, and weighting so that we can append logfree next
pop_module(); pop_module();

% Stop fitting the PZ wavelet.
STACK{end-2}{1}.fit_fields = {};

end
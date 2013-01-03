function linear_fit_spn(cellid, training_set)
% Fits a linear model directly on the output of envelopes provided by BAPHY

global NARF_PATH STACK XXX;
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

XXX = {};
XXX{1} = [];
XXX{1}.cellid = cellid;
XXX{1}.training_set = training_set;
XXX{1}.test_set = {};

raster_fs = 200;
filter_length = 20;
n_channels = 2;

STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy.mdl(...
                struct('raw_resp_fs', raster_fs, ...
                       'raw_stim_fs', raster_fs,...
                       'stimulus_format','envelope', ...
                       'stimulus_channel_count', n_channels));
STACK{2} = mdls.fir_filter.mdl(struct('num_dims', n_channels, ...
                                      'num_coefs', filter_length));

recalc_xxx(1);

STACK{2}.fit_fields = {'coefs'};

fit_with_lsqcurvefit();

STACK{3} = mdls.correlation;
STACK{4} = mdls.mean_squared_error;

recalc_xxx(2);  % Recompute now from the FIR filter onward


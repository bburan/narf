function linear_fit_spn_volterra(cellid, training_set)
% Fits a volterra model with a polynomial output linearity.

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
STACK{2} = mdls.concat_second_order_terms;                   
STACK{3} = mdls.normalize_channels;                  
STACK{4} = mdls.fir_filter.mdl(...
                struct('num_dims', n_channels + nchoosek(n_channels, 2), ...
                       'num_coefs', filter_length));
STACK{5} = mdls.nonlinearity.mdl(struct('phi', [1 1 1 1], ...
                                        'nlfn', @polyval));

recalc_xxx(1);

STACK{4}.fit_fields = {'coefs'};
STACK{5}.fit_fields = {'phi'};

fit_with_lsqcurvefit();

STACK{6} = mdls.mean_squared_error;
STACK{7} = mdls.correlation;

recalc_xxx(2); 

function linear_fit_spn_exp_v2(cellid, training_set)
% Fits a linear model with an exponential output
% Trains in three batches.

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
STACK{2} = mdls.normalize_channels;                   
STACK{3} = mdls.fir_filter.mdl(struct('num_dims', n_channels, ...
                                      'num_coefs', filter_length));

% Fit the coefficients first
recalc_xxx(1);
STACK{3}.fit_fields = {'coefs'};
fit_with_lsqcurvefit();                
                                 
% Then add the nonlinearity
STACK{4} = mdls.nonlinearity.mdl(struct('phi', [1 0], ...
                                        'nlfn', @exponential)); 
recalc_xxx(3); 

% Then fit the exponential
STACK{3}.fit_fields = {};
STACK{4}.fit_fields = {'phi'};
fit_with_lsqcurvefit();

% Finally, fit all parameters
STACK{3}.fit_fields = {'coefs'};
STACK{4}.fit_fields = {'phi'};
fit_with_lsqcurvefit();

% Now the reporting
STACK{5} = mdls.correlation;

recalc_xxx(2); 


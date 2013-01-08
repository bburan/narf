function linear_fit_spn_exp_sigmoid(cellid, training_set)
% Fits a linear model with a sigmoid output nonlinearity

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
            
STACK{3}.fit_fields = {'coefs'};
fit_with_lsqcurvefit();             


STACK{4} = mdls.nonlinearity.mdl(struct('phi', [0.05 0.05 0.1 0], ...
                                        'nlfn', @sigmoidal));    


STACK{4}.fit_fields = {'phi'};
fit_with_lsqcurvefit();

STACK{5} = mdls.nonlinearity.mdl(struct('phi', [1 0], ...
                                        'nlfn', @exponential));    
STACK{5}.fit_fields = {'phi'};
fit_with_lsqcurvefit();



% Now the reporting
STACK{6} = mdls.correlation;

recalc_xxx(3); 

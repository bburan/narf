function linear_fit_spn_stephen(cellid, training_set)
% Fits a linear model directly on the output of envelopes provided by BAPHY

fprintf('Fitting linear model to SPN: %s\n', cellid);

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

                                  
recalc_xxx(1);

m = [];
m.num_dims = n_channels;
m.altcore = 'cdcore';
m.maxlag = filter_length - 1;
m.resampcount = filter_length - 1;  
m.sfscount = 10;
m.sfsstep = 3;
m.rasterfs = raster_fs;
m.coefs = zeros(n_channels, filter_length);
m.stimfield = 'stim'; % Names of the fields to use as stimulus and response
m.respfield = 'resp';

STACK{3}.coefs = do_stephen_fit(m, XXX{3});
ss = size(STACK{3}.coefs);
STACK{3}.num_dims = ss(1);
STACK{3}.num_coefs = ss(2);

% Add a nonlinear stage
STACK{4} = mdls.nonlinearity.mdl(struct('phi', [1 1 1 1], ...
                                        'nlfn', @polyval));                       
STACK{4}.fit_fields = {'phi'};

recalc_xxx(2);     
fit_with_lsqcurvefit();

STACK{5} = mdls.correlation;

recalc_xxx(2);


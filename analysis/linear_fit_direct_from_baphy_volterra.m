% Fits a linear model directly on the output of envelopes provided by BAPHY

narf_set_path;

global NARF_PATH STACK XXX;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

% XXX defines the initial data
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por024a-a1';
XXX{1}.training_set = {'por024a19_p_SPN'};
XXX{1}.test_set = {};

% Define the model structure, using defaults for unspecified fields.
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
STACK{3} = mdls.fir_filter.mdl(...
                struct('num_dims', n_channels + nchoosek(n_channels, 2), ...
                       'num_coefs', filter_length));
                  
STACK{4} = mdls.nonlinearity.mdl(struct('nlfn', @exponential, ...
                                        'phi', [1 1]));
STACK{5} = mdls.correlation;
STACK{6} = mdls.mean_squared_error;

% Compute the whole XXX structure once
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

% Run Stephen's fitting routine and update stack from that point on
STACK{3}.coefs = do_stephen_fit(m, XXX{3});
ss = size(STACK{3}.coefs);
STACK{3}.num_dims = ss(1);
STACK{3}.num_coefs = ss(2);
recalc_xxx(3); 

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 900]);
narf_modelpane(pf, mdls); 

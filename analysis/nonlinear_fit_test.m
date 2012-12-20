% Test of how to use lsqcurvefit to fit a polynomial output nonlinearity

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

% TODO: Do a linear fit here if you like
tmp = load('volterra_coefs.mat');
STACK{3} = mdls.fir_filter.mdl(...
                struct('num_dims', n_channels + nchoosek(n_channels, 2), ...
                       'num_coefs', filter_length, ...
                       'coefs', tmp.volterra_coefs));
STACK{4} = mdls.nonlinearity.mdl(struct('phi', [0 0 0], ...
                                        'nlfn', @polyval)); 


% Compute the whole XXX structure once
recalc_xxx(1);

% Fit nonlinear coefficients
STACK{4}.fit_fields = {'phi'};
phi_init = pack_fittables(STACK);
% options = optimset('Display','final','LargeScale','off');
%options = optimset('DiffMinChange',1e1);
phi_best = lsqcurvefit(@my_mse, phi_init, 4, XXX{4}.dat.por024a19_p_SPN.respavg(:), [], []);
unpack_fittables(phi_best);

% Now compute the MSE and Correlations
STACK{5} = mdls.mean_squared_error;
STACK{6} = mdls.correlation;

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 900]);
narf_modelpane(pf, mdls); 

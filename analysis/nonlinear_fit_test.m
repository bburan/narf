% Test of how to use lsqcurvefit to fit a polynomial output nonlinearity

narf_set_path;

global NARF_PATH STACK XXX;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

% XXX defines the initial data
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por026b-b2'; %'por024a-a1';
XXX{1}.training_set = {'por026b09_p_SPN'}; %{'por024a19_p_SPN'};
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
STACK{3} = mdls.normalize_channels;                  

% Do an initial linear fit here if you like...or load some initial values
tmp = load('volterra_coefs.mat');
STACK{4} = mdls.fir_filter.mdl(...
                struct('num_dims', n_channels + nchoosek(n_channels, 2), ...
                       'num_coefs', filter_length, ...
                       'coefs', tmp.volterra_coefs));
% Polynomial
% STACK{5} = mdls.nonlinearity.mdl(struct('phi', [0 0 0], ...
%                                         'nlfn', @polyval)); 

% Sigmoidal
% STACK{5} = mdls.nonlinearity.mdl(struct('phi', [1 pi/2 0.1], ...
%                                         'nlfn', @sigmoidal)); 

% Exponential
STACK{5} = mdls.nonlinearity.mdl(struct('phi', [1 0], ...
                                        'nlfn', @exponential)); 

% Dumb thresholding and scaling
% STACK{5} = mdls.nonlinearity.mdl(...
%     struct('phi', [1 0], ...
%            'nlfn', @(phi, z) (z + phi(2)) .* phi(1) .* heaviside(z + phi(2)))); 

% Compute the whole XXX structure once to update it
recalc_xxx(1);

% Fit nonlinear coefficients and linear coefficients at same time
STACK{4}.fit_fields = {'coefs'};
STACK{5}.fit_fields = {'phi'};

fit_with_lsqcurvefit();

% Now compute the MSE and Correlations
STACK{6} = mdls.mean_squared_error;
STACK{7} = mdls.correlation;

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 900]);
narf_modelpane(pf, mdls); 

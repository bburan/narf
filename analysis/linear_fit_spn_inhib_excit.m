function linear_fit_spn_inhib_excit(cellid, training_set)
% Fits an inhibition_excitation model with polynomial nonlinearities:
%
%             /---> FIR ---> RAMP FN ----> POLYNOMIAL --\
%  SIGNAL ---<                                          (+)---> OUTPUT
%             \---> FIR ---> RAMP FN ----> POLYNOMIAL --/
%

global NARF_PATH STACK XXX;
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

XXX = {};
XXX{1} = [];
XXX{1}.cellid = cellid;
XXX{1}.training_set = training_set;
XXX{1}.test_set = {};

raster_fs = 200;
filter_length = 10;
n_channels = 2;

STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy.mdl(...
                struct('raw_resp_fs', raster_fs, ...
                       'raw_stim_fs', raster_fs,...
                       'stimulus_format','envelope', ...
                       'stimulus_channel_count', n_channels));
STACK{2} = mdls.normalize_channels;

% The two filters
STACK{3} = mdls.fir_filter.mdl(...
                struct('num_dims', n_channels, ...
                       'num_coefs', filter_length, ...
                       'output', 'inhib'));
STACK{4} = mdls.fir_filter.mdl(...
                struct('num_dims', n_channels, ...
                       'num_coefs', filter_length, ...
        		       'output', 'excit'));

% Try to make convergence occur with positive coefficients
STACK{3}.coefs = 0.01*ones(size(STACK{3}.coefs));
STACK{4}.coefs = 0.01*ones(size(STACK{3}.coefs));                   
                   

% % Guarantee the outputs will be only one-sided
% STACK{5} = mdls.nonlinearity.mdl(struct('input_stim', 'inhib', ...
% 					'output', 'inhib', ...
% 					'phi', [-1 0], ...
%                     'nlfn', @(phi, z) (z + phi(2)) .* phi(1) .* heaviside(z + phi(2))));
% 
% STACK{6} = mdls.nonlinearity.mdl(struct('input_stim', 'excit', ...
% 					'output', 'excit', ...
% 					'phi', [1 0], ...
%                     'nlfn', @(phi, z) (z + phi(2)) .* phi(1) .* heaviside(z + phi(2))));

% Put sigmoidal curves on the filter outputs
STACK{5} = mdls.nonlinearity.mdl(struct('input_stim', 'inhib', ...
                                        'output', 'inhib', ...
                                        'phi', [0.05 0.05 -0.1 0], ...
                                        'nlfn', @sigmoidal));
                                    
STACK{6} = mdls.nonlinearity.mdl(struct('input_stim', 'excit', ...
                                        'output', 'excit', ...
                                        'phi', [0.05 0.05 0.1 0], ...
                                        'nlfn', @sigmoidal));

STACK{7} = mdls.sum_fields.mdl(struct('inputs', {{'inhib', 'excit'}}, ...
                                      'output', 'stim'));			   

recalc_xxx(1);

STACK{3}.fit_fields = {'coefs'};
STACK{4}.fit_fields = {'coefs'};
STACK{5}.fit_fields = {'phi'};
STACK{6}.fit_fields = {'phi'};

fit_with_lsqcurvefit();

STACK{8} = mdls.correlation;

recalc_xxx(3); 

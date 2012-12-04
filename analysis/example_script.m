% A script to demonstrate how to use NARF and 'narf_modelpane.m'
%
% Essentially, this script sets up a model using a gammatone filter bank
% prefilter, a simple decimation using the mean value and square root, and
% an FIR filter which accepts the output of every gammatone filter. When
% you feed it some data, boosting is used to estimate the parameters of
% a FIR filter. The FIR filter coefficients should then define an STRF.
%
% Note that this is a very primitive first attempt, and should not be used
% for any real scientific use, mostly because the boosting routine is
% really really simple and naive.

narf_set_path;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH 'modules/']);

global STACK XXX;

% XXX defines the initial data
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por012c-b1';
XXX{1}.training_set = {'por012c02_p_TOR'};
XXX{1}.test_set = {};

% Define the model structure
STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy;
STACK{2} = mdls.gammatone_filter_bank.mdl(struct('num_channels', 10));
STACK{3} = mdls.downsample_with_fn;
STACK{4} = mdls.fir_filter.mdl(struct('num_dims', 10, ...
                                      'fit_fields', {'coefs'}));

% TODO: Trigger a computation to update the stack
                                  
% Define specific parameter values to try (in a non-empty cell array) 
% Or search through other parameters freely (in anything else)

%= cell(size(STACK));
%{3, 'preconf_fn', {@min, @mean, @max, @sum}};
%{4, 'postconf_fn', {@(x) x, @sqrt, @(x) log(x+10^-3)}};

% Define the free parameters which optimization may use for anything


% Run the optimization for each of the iterated dimensions


% Save the stack
% save_stack(STACK, 'example.m')

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 600]);

narf_modelpane(pf, mods); 

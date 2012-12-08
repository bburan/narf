% A script to demonstrate how to compare models with discrete values for
% certain parameters.

% Define specific parameter values to try (in a non-empty cell array) 
% Or search through other parameters freely (in anything else)
%{3, 'preconf_fn', {@min, @mean, @max}};
%{4, 'postconf_fn', {@(x) x, @sqrt, @(x) log(x+10^-3)}};

% Run the optimization for each of the iterated dimensions
narf_set_path;

global NARF_PATH STACK XXX;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH 'modules/']);

% XXX defines the initial data
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por022b-a1';
XXX{1}.training_set = {'por022b02_p_TOR'};
XXX{1}.test_set = {};

% Define the model structure, using defaults for unspecified fields.
STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy;
STACK{2} = mdls.gammatone_filter_bank.mdl(struct('num_channels', 10, ...
                                                 'bank_min_freq', 250, ...
                                                 'bank_max_freq', 8000));
STACK{3} = mdls.downsample_with_fn;
STACK{4} = mdls.fir_filter.mdl(struct('num_dims', 10, ...
                                      'num_coefs', 20));
                                  
% NOTE: Must use double braces {{ }} because of struct()'s default behavior
% of unfortunately stripping away the first layer of cell references.

% Recompute the entire stack
recalc_stack(1);

m = [];
m.num_dims = 2;
m.altcore = 'cdcore';
m.maxlag = 19;
m.resampcount = 19;  % Why is the total number of things this + 1?
m.sfscount = 10;
m.sfsstep = 3;
m.rasterfs = 200;  % TODO: REMOVE ME
m.coefs = zeros(m.num_dims, 10);

% Run Stephen's fitting routine at depth 4
STACK{4}.coefs = do_stephen_fit(m, XXX{4});
ss = size(STACK{4}.coefs);
STACK{4}.num_dims = ss(1);
STACK{4}.num_coefs = ss(2);
% Recompute now just the last little bit
recalc_stack(4);

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 600]);
narf_modelpane(pf, mdls); 

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
STACK{2} = mdls.gammatone_filter_bank.mdl(struct('num_channels', 2, ...
                                                 'bank_min_freq', 250, ...
                                                 'bank_max_freq', 8000));
STACK{3} = mdls.downsample_with_fn;
STACK{4} = mdls.fir_filter.mdl(struct('num_dims', 2, ...
                                      'fit_fields', {{'coefs'}}));
% NOTE: Must use double braces {{ }} because of struct()'s default behavior
% of unfortunately stripping away the first layer of cell references.

% Run test optimization routine (stupid boost that stops on iteration #)
n_iters = 5;
x_0 = pack_fittables(STACK);
if isempty(x_0)
    log_msg('No parameters were selected to be fit.');
    return;
end
[x_bst, s_bst] = boosting(x_0', @correlation_of_downsampled_signals, ...
                         @(n,x,s)(n > n_iters), 1.0);

unpack_fittables(x_bst);
 
% coefs = STACK{end}.coefs;
% save('real_good_coefs.mat', 'coefs');

% Save the model stack somewhere
% save_model(STACK, 'example.m')

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 600]);
narf_modelpane(pf, mdls); 

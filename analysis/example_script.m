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

global NARF_PATH STACK XXX;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH 'modules/']);

% XXX defines the initial data
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por012c-b1';
XXX{1}.training_set = {'por012c02_p_TOR'};
XXX{1}.test_set = {};

% Define the model structure
STACK = {};
STACK{1} = mdls.load_stim_resps_from_baphy;
STACK{2} = mdls.gammatone_filter_bank.mdl(struct('num_channels', 2));
STACK{3} = mdls.downsample_with_fn;
STACK{4} = mdls.fir_filter.mdl(struct('num_dims', 2, ...
                                      'fit_fields', {{'coefs'}}));
% NOTE: Must use double braces {{ }} because of struct()'s default behavior
% which strips away the first layer of cell references. 
                                  
% Run test optimization routine (stupid boost that stops on iteration #)
n_iters = 10;
x_0 = pack_fittables(STACK);
if isempty(x_0)
    log_msg('No parameters were selected to be fit.');
    return;
end
[x_bst, s_bst] = boosting(x_0', @correlation_of_downsampled_signals, @(n,x,s)(n > n_iters), 1.0);
unpack_fittables(x_bst);
update_tables_and_plots();

% Save the stack
% save_stack(STACK, 'example.m')

% Finally, display the GUI for easy tweaking and viewing of the best result
pf = figure('Menubar','figure', 'Resize','off', ...
            'Units','pixels', 'Position', [20 50 1300 600]);
narf_modelpane(pf, mdls); 

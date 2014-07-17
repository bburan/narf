function lognzn()
% With base-n log and a zero-below threshold

global MODULES STACK;


append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

nl_log_nozero = @(p,z) nl_log([p(1) 0 p(2)], z);

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 0], ...
                                              'nlfn', nl_log_nozero, ...
                                              'fit_fields', {{'phi'}}, ...
                                              'fit_constraints', {{struct( ...
                                                       'var', 'phi', ...
                                                       'lower', [-50 -50], ...
                                                       'upper', [50 50])}} ...
                                               )));
STACK{end}{1}.auto_plot=[];

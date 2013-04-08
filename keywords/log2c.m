function log2c()
% Adds a log2 zero baseline channel as a NEW channel without overwriting
% the original signal.

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], ...
                                              'nlfn', @nl_log, ...
                                              'output', 'stimlog2c')));
                                                                                    
append_module(MODULES.concatenate_channels.mdl(struct('input1', 'stim', ...
                                                      'input2', 'stimlog2c')));
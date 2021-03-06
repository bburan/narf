% 100 hz sampling, sqrt of response vectors
function env100rr()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 100, ...
                                  'raw_stim_fs', 100,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'envelope'))); 

append_module(MODULES.nonlinearity.mdl(struct('phi', [2], ...
                                              'nlfn', @nl_root,...
                                              'input_stim','respavg',...
                                              'output','respavg')));

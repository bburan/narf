function env400()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 400, ...
                                  'raw_stim_fs', 400,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'envelope'))); 
                              

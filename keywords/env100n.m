function env100n()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 100, ...
                                  'raw_stim_fs', 100, ...
                                  'include_prestim', false,...
                                  'stimulus_format', 'envelope'))); 
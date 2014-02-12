function gtSNS30()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 100, ...
                                  'raw_stim_fs', 100,...
                                  'stimulus_channel_count', 30, ...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'gamma'))); 

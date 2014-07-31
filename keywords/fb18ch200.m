function fb18ch200()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 200, ...
                                  'raw_stim_fs', 200,...
                                  'stimulus_channel_count', 18, ...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'ozgf'))); 

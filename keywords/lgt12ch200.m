function lgt12ch200()

global MODULES;
global ANESTHETIZED_MOUSE

ANESTHETIZED_MOUSE=1;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 200, ...
                                  'raw_stim_fs', 200,...
                                  'stimulus_channel_count', 12, ...
                                  'include_prestim', 1, ...
                                  'stimulus_format', 'gamma264', ...
                                  'response_format', 'lfp+pupil'))); 

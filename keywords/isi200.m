function isi200()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 10000, ...
                                  'raw_stim_fs', 200,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'envelope'))); 

append_module(MODULES.inter_spike_intervals); 

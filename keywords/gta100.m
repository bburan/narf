function gta100()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 100, ...
                                  'raw_stim_fs', 100000,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'wav'))); 
                              
append_module(MODULES.gammatone_filter_bank.mdl(...
                           struct('bank_min_fs', 500, ...
                                  'bank_max_fs', 5000, ...
                                  'raw_stim_fs', 100000,...
                                  'use_env', false, ...
                                  'align_phase', false, ...
                                  'num_channels', 10))); 
                              
append_module(MODULES.downsample_signal.mdl(...
                    struct('input_freq', 100000, ...
                           'input', 'stim', ...
                           'input_time', 'stim_time', ...
                           'output_freq', 100, ...
                           'output', 'stim', ...
                           'output_time', 'stim_time'))); 

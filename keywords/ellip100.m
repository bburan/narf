function ellip100()

global MODULES;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                                      struct('raw_resp_fs', 200, ...
                                             'raw_stim_fs', 100000)));
                                         
append_module(MODULES.elliptic_bandpass_filter_bank);


append_module(MODULES.downsample_with_fn.mdl(...
                                    struct('downsampled_freq', 200, ...
                                           'conv_fn', @mean, ...
                                           'postconv_fn', @abs)));
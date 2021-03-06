function fft128ch200()

global MODULES STACK XXX META;

channels = 128;
load_fs = 50000;
samp_fs = 200;

append_module(MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', samp_fs, ...
                                  'raw_stim_fs', load_fs,...
                                  'include_prestim', true, ...
                                  'stimulus_format', 'wav'))); 

append_module(MODULES.fourier_transform.mdl(...
                struct('n_output_chans', channels, ...
                       'min_freq', 200, ...
                       'input_freq', load_fs, ...                                              
                       'output', 'stim', ...
                       'output_freq', samp_fs)));  

normp();
%wc01();
%nmse();
%fit05();
%fit_lsq();

%fitell();
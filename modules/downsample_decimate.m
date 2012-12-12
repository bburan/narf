function m = downsample_decimate(args)
% Wraps the behavior of MATLAB's decimate() function.
%
% Module fields that must ALWAYS be defined
m = [];
m.mdl = @downsample_decimate;
m.name = 'downsample_decimate';
m.fn = @do_downsampling;
m.pretty_name = 'Use MATLAB''s decimate()';
m.editable_fields = {'downsampled_freq', 'pre_ds_fn', 'post_ds_fn', ...
                     'input', 'input_time', 'output', 'output_time'};
m.isready_pred = @downsampler_isready;

% Module fields that are specific to THIS MODULE
m.downsampled_freq = 200;
m.pre_ds_fn = @abs;
m.post_ds_fn = @sqrt;
m.input = 'pp_stim';
m.input_time = 'raw_stim_time';
m.output = 'ds_stim';
m.output_time = 'ds_stim_time';
m.output_fs = 'ds_stim_fs';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_downsampled_stimulus;
m.plot_fns{1}.pretty_name = 'Downsampled Stimulus vs Time';

function x = do_downsampling(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    scale = floor(baphy_mod.raw_stim_fs / mdl.downsampled_freq);

    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, N, F] = size(x.dat.(sf).(mdl.input));
        x.dat.(sf).(mdl.output) = zeros(S,ceil(N/scale),F);
        for s = 1:S
            for f = 1:F
                x.dat.(sf).(mdl.output)(s,:,f) = m.post_ds_fn(...
                    decimate(m.pre_ds_fn(x.dat.(sf).(mdl.input)(s,:,f)), scale));      
            end    
        end
        x.dat.(sf).(mdl.output_time) = ...
            linspace(1/mdl.downsampled_freq, ...
                     x.dat.(sf).(mdl.input_time)(end), ...
                     length(x.dat.(sf).(mdl.output)));
        x.dat.(sf).(mdl.output_fs) = mdl.downsampled_freq; %% TODO: Remove SPOT violation
  
    end
                     
end

end
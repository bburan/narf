function m = downsample_matlab(args)
% Wraps the behavior of MATLAB's downsample() function.
%
% Module fields that must ALWAYS be defined
m = [];
m.mdl = @downsample_matlab;
m.name = 'downsample_matlab';
m.fn = @do_downsampling;
m.pretty_name = 'Use MATLAB''s downsample()';
m.editable_fields = {'downsampled_freq', 'pre_ds_fn', 'post_ds_fn', ...
                     'input', 'input_time', 'output', 'output_time'};
m.isready_pred = @isready_always;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_channel_vs_time(stack, xxx, m.output_time, m.output);
m.plot_fns{1}.pretty_name = 'Downsampled Channel vs Time';

% Module fields that are specific to THIS MODULE
m.downsampled_freq = 200;
m.pre_ds_fn = @abs;
m.post_ds_fn = @sqrt;
m.input = 'stim';
m.input_time = 'stim_time';
m.output = 'stim';
m.output_time = 'stim_time';
m.output_fs = 'stim_fs';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_downsampling(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    scale = floor(baphy_mod.raw_stim_fs / mdl.downsampled_freq);

    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, N, C, F] = size(x.dat.(sf).(mdl.input));
        x.dat.(sf).(mdl.output) = zeros(S,ceil(N/scale),F);
        for s = 1:S
            for c = 1:C
                for f = 1:F
                    x.dat.(sf).(mdl.output)(:,s,c,f) = m.post_ds_fn(...
                        downsample(m.pre_ds_fn(x.dat.(sf).(mdl.input)(:,s,c,f)), scale));      
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

end
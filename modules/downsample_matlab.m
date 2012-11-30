function m = downsample_matlab(args)
% Wraps the behavior of MATLAB's downsample() function.
%
% Module fields that must ALWAYS be defined
m = [];
m.mdl = @downsample_matlab;
m.name = 'downsample_matlab';
m.fn = @do_downsampling;
m.pretty_name = 'Use MATLAB''s downsample()';
m.editable_fields = {'downsampled_freq', 'pre_ds_fn', 'post_ds_fn'};
m.isready_pred = @downsampler_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_downsampled_stimulus;
m.plot_fns{1}.pretty_name = 'Downsampled Stimulus vs Time';

% Module fields that are specific to THIS MODULE
m.downsampled_freq = 200;
m.pre_ds_fn = @abs;
m.post_ds_fn = @sqrt;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_downsampling(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    scale = ceil(baphy_mod.raw_stim_fs / mdl.downsampled_freq);

    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, N, F] = size(x.dat.(sf).pp_stim);
        x.dat.(sf).ds_stim = zeros(S,N/scale,F);
        for s = 1:S
            for f = 1:F
                x.dat.(sf).ds_stim(s,:,f) = m.post_ds_fn(...
                    downsample(m.pre_ds_fn(x.dat.(sf).pp_stim(s,:,f)), scale));      
            end    
        end
        x.dat.(sf).ds_stim_time = ...
                    linspace(1/mdl.downsampled_freq, ...
                             x.dat.(sf).raw_stim_time(end), ...
                             length(x.dat.(sf).ds_stim));
                         
        x.dat.(sf).ds_stim_fs = mdl.downsampled_freq;
    end

end

end
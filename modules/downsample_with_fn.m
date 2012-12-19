function m = downsample_with_fn(args)
% A function to create a (possibly nonlinear) downsampling module.
%
% Very simple algorithm to downsample a signal:
%    Apply PRECONV_FN to every value of the input signal.
%    Take N of those values at a time, where N = ceil(raw_freq / downsampled_freq)
%    Apply those N values to function CONV_FN. (Do the convolution, baby!)
%    Apply POSTCONV_FN to each of those values. 
%
% The only PRECONV_FN function I can think is useful is @abs, but maybe
% somebody will want to use another function someday.
%
% Good CONV_FN functions include @max, @min, @mean, and @sum, but you may
% prefer to define your own.
%
% Good POSTCONV_FN functions are @sqrt, and maybe @log if your data never
% goes to zero. In such a case @(x) log(x+10^-6) may be nice, substituting
% the 10^-6 with your noise floor intensity.
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how MODULES are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @downsample_with_fn;
m.name = 'downsample_with_fn';
m.fn = @do_downsampling;
m.pretty_name = 'Downsample with Arbitrary Function';
m.editable_fields = {'downsampled_freq', ...
                     'preconv_fn', 'conv_fn', 'postconv_fn', ...
                     'input', 'input_time', 'output', 'output_time'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.downsampled_freq = 200;
m.preconv_fn = @abs;
m.conv_fn = @mean;
m.postconv_fn = @sqrt;
m.input = 'stim';
m.input_time = 'stim_time';
m.output = 'stim';
m.output_time = 'stim_time';
m.output_fs = 'stim_fs';

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
        fprintf('Downsampling %s\n', sf);
        x.dat.(sf).(mdl.output) = ...
            mdl.postconv_fn(conv_fn(mdl.preconv_fn(x.dat.(sf).(mdl.input)), ...
                                    1, mdl.conv_fn, scale, 0));
        x.dat.(sf).(mdl.output_time) = ...
            linspace(1/mdl.downsampled_freq, ...
                     x.dat.(sf).(mdl.input_time)(end), ...
                     length(x.dat.(sf).(mdl.output)));
        x.dat.(sf).(mdl.output_fs) = mdl.downsampled_freq; %% TODO: Remove SPOT violation
    end
end

end
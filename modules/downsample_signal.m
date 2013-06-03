function m = downsample_signal(args)
% Downsamples a signal using matlab's downsample() function
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how MODULES are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @downsample_signal;
m.name = 'downsample_signal';
m.fn = @do_downsampling;
m.pretty_name = 'Downsample Signal';
m.editable_fields = {'input', 'input_freq', 'input_time', ...
                     'output', 'output_freq', 'output_time', 'downsampler'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'respavg';
m.input_freq = 10000;
m.input_time = 'resp_time';
m.output = 'respavg';
m.output_freq = 200;
m.output_time = 'respavg_time';
m.downsampler = @downsample; % Also try @decimate, or @conv_fn, or @resample

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_downsampled;
m.plot_fns{1}.pretty_name = 'Output vs Time';

function x = do_downsampling(mdl, x, stack, xxx)
    
    scale = mdl.input_freq / mdl.output_freq;
    
    if ~(floor(scale) == scale)
        error('Input frequency must be an integer multiple of output frequency.');
    end

    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        out = zeros(ceil(T/scale),S,C);
        for s = 1:S
            for c = 1:C                
                out(:,s,c) = mdl.downsampler(x.dat.(sf).(mdl.input)(:,s,c), scale);
            end
        end
        x.dat.(sf).(mdl.output) = out;                 
        
        x.dat.(sf).(mdl.output_time) = ...
                linspace(1/mdl.output_freq, ...
                     x.dat.(sf).(mdl.input_time)(end), ...
                     length(x.dat.(sf).(mdl.output)))';        
    end
end

function do_plot_downsampled(sel, stack, xxx)
    % [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    mdls = stack{end};
    xouts = {xxx{end}};
    
    hold on;
    do_plot(xouts, mdls{1}.output_time, mdls{1}.output, ...
            sel, 'Time [s]', 'Respavg [-]');
    hold off;     
end

end
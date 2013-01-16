function m = passthru(args)
% A module that does not do anything, just passes the signals on thru.
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @passthru;
m.name = 'passthru';
m.fn = @do_passthru;
m.pretty_name = 'Pass-Thru';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input  = 'stim';
m.time   = 'stim_time';
m.output = 'stim'; 

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_channel_vs_time(stack, xxx, m.time, m.input);
m.plot_fns{1}.pretty_name = 'Input Channel vs Time';
m.plot_gui_create_fn = @create_chan_selector_gui;

function x = do_passthru(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    sfs = fieldnames(x.dat);
    for ii = 1:length(x.dat)
        sf = sfs{ii};
        x.dat.(sf).(mdl.output) = x.dat.(sf).(mdl.input);
    end
end
end
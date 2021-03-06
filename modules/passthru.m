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

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

% Optional fields
m.plot_fns = {};

function x = do_passthru(mdl, x)
    if ~strcmp(mdl.input, mdl.output)
        fields = fieldnames(x.dat);
        for ii = 1:length(fields)
            sf = fields{ii};
            x.dat.(sf).(mdl.output) = x.dat.(sf).(mdl.input);
        end
    end
end

end

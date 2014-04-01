function m = hilbert_transform(args)
% Hilbert Transform. Takes the signal envelope.
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @hilbert_transform;
m.name = 'hilbert_transform';
m.fn = @do_hilbert_transform;
m.pretty_name = 'Hilbert Transform';
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
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Output Channel';

function x = do_hilbert_transform(mdl, x, stack, xxx)
    sfs = fieldnames(x.dat);    
    for ii = 1:length(x.dat)
        sf = sfs{ii};        
        x.dat.(sf).(mdl.output) = abs(hilbert(x.dat.(sf).(mdl.input)));
    end
end

end
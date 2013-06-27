function m = concatenate_channels(args)
% Concatenates two signals together

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @concatenate_channels;
m.name = 'concatenate_channels';
m.fn = @do_concatenate_channels;
m.pretty_name = 'Concatenate Channels';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim'; 
m.input2 = 'stim2';
m.output = 'stim';

% Optional fields
m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Channels (All)';
m.plot_fns{2}.fn = @do_plot_single_default_output;
m.plot_fns{2}.pretty_name = 'Channel (Single)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_concatenate_channels(mdl, x, stack, xxx)    

    for sf = fieldnames(x.dat)', sf=sf{1};           
        [T1, S1, C1] = size(x.dat.(sf).(mdl.input1));
        [T2, S2, C2] = size(x.dat.(sf).(mdl.input1));
        
        if (T1 ~= T2 || S1 ~= S2)
            error('Concatenated channels must have the same length in time and stimulus dimensions!');
        end

        out = zeros(T1,S1, C1+C2);
        
        for ii = 1:(C1+C2) 
            if (ii <= C1)
                out(:,:,ii) = x.dat.(sf).(mdl.input1)(:,:,ii);
            else
                out(:,:,ii) = x.dat.(sf).(mdl.input1)(:,:,ii-C1);
            end
        end
        
        x.dat.(sf).(mdl.output) = out;
    end
end

end
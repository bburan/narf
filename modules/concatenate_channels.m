function m = concatenate_channels(args)
% Concatenates multiple signals together

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @concatenate_channels;
m.name = 'concatenate_channels';
m.fn = @do_concatenate_channels;
m.pretty_name = 'Concatenate Channels';
m.editable_fields = {'inputs', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.inputs = {'stim1', 'stim2'}; 
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

        % Allocate space
        Ts=[]; Ss=[]; Cs=[];        
        for cc = 1:length(mdl.inputs)
            [T, S, C] = size(x.dat.(sf).(mdl.inputs{cc}));            
            Ts(end+1) = T;
            Ss(end+1) = S;
            Cs(end+1) = C;
        end
        if all(Ts(1) == Ts) && all(Ss(1) == Ss)
            out = zeros(max(Ts),max(Ss), sum(Cs));
        else
            fprintf('Ts = %d', Ts);
            fprintf('Ss = %d', Ss);            
            error('Error in sizes of concatenated channels\n');
        end
    
        % Concatenate all channels
        c_idx = 1;
        for cc = 1:length(mdl.inputs)
            [T, S, C] = size(x.dat.(sf).(mdl.inputs{cc}));            
            for ii = 1:C 
                out(:,:,c_idx) = x.dat.(sf).(mdl.inputs{cc})(:,:,ii);
                c_idx = c_idx + 1;
            end
        end
        x.dat.(sf).(mdl.output) = out;
    end
end

end
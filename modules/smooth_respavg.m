function m = smooth_respavg(args)
% Smooths the respavg signal ONLY in the training set. 
% Uses an arbitrarily defined kernel.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @smooth_respavg;
m.name = 'smooth_respavg';
m.fn = @do_smooth_respavg;
m.pretty_name = 'Smooth Respavg';
m.editable_fields = {'input', 'time', 'output', 'window'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'respavg'; 
m.time  = 'resp_time';
m.output = 'respavg';
m.window = [1 4 1];

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_respavg;
m.plot_fns{1}.pretty_name = 'Training Set Respavg';
m.plot_gui_create_fn = @create_chan_selector_gui;

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_smooth_respavg(mdl, x, stack, xxx)

    % Build the gaussian smoothing filter   
    % gf = gausswin(mdl.window_size);
    gf = mdl.window;
    gf = gf / sum(gf);    
    
    % Smooth ONLY the training set RESPAVG signals. 
    for ii = 1:length(x.training_set)
        sf = x.training_set{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        out = zeros([T, S, C]);
        
        for c = 1:C,
            for s = 1:S
                out(:,s,c) = conv(x.dat.(sf).(mdl.input)(:,s,c), gf, 'same');
            end
        end
        x.dat.(sf).(mdl.output) = out;
    end
end

function do_plot_respavg(sel, stack, xxx)
    mdl = stack{end}{1};
    xold = xxx{end-1};
    xnew = xxx{end};
    
    %[sf, stim_idx, baphy_chan_idx] = get_baphy_plot_controls(stack);
    %chan_idx = popup2num(mdl.plot_gui.selected_chan_popup);
    
    plot(xold.dat.(sel.stimfile).(mdl.time), ...
         xold.dat.(sel.stimfile).(mdl.input)(:, sel.stim_idx, sel.chan_idx), 'k-', ...
         xold.dat.(sel.stimfile).(mdl.time), ...
         xnew.dat.(sel.stimfile).(mdl.output)(:, sel.stim_idx, sel.chan_idx), 'r-');
     
end

end
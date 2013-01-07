function hs = create_chan_selector_gui(parent_handle, stack, xxx)
    % Please run AFTER execution so that output signal exists
    pos = get(parent_handle, 'Position');
    w = pos(3) - 10;
    h = pos(4) - 10;
    hs = [];
    
    m = stack{end};
    mod_idx = length(stack);
    x = xxx{end};
    
    % Create a channel selector
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'Chan:', ...
        'Units', 'pixels', 'Position', [5 (h-25) 50 25]);
    hs.selected_chan_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [45 (h-25) w-50 25], ...
        'Callback', @(a,b,c) selected_chan_popup_callback());

    function update_channel_popup(stack)
        % Find the baphy module's selected stimfile
        sf = popup2str(find_module_gui_control(stack, 'selected_stimfile_popup'));
        
        if isfield(x.dat, sf)
            [T, S, C] = size(x.dat.(sf).(m.output));
            d = {};
            for ii = 1:C
                d{ii} = sprintf('%d',ii);
            end
            set(hs.selected_chan_popup, 'String', char(d));
            set(hs.selected_chan_popup, 'Value', 1);
        else
            error('Selected stimulus file not found: %s', sf);
        end
    end
    
    function selected_chan_popup_callback()
        % Call the plot function again via the plot_popup 
        hgfeval(get(m.gh.plot_popup, 'Callback'), mod_idx, []);
        drawnow;
    end

    update_channel_popup(stack);

end
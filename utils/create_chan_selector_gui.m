function hs = create_chan_selector_gui(parent_handle, stack, xxx, chanfield)
% hs = create_chan_selector_gui(parent_handle, stack, xxx, chanfield)
%
% Creates a UICONTROL popup menu in parent_handle under the fieldname
% 'selected_chan_popup'. This can then be used by default methods like
% do_plot_channel_vs_time() so that the amount of boilerplate
% copy-and-paste code for a module can be reduced since they share the same
% implementation. 
%
% ARGUMENTS:
%   parent_handle       The panel handle in which the popup will be placed
%
%   stack               Local copy of STACK
%   xxx                 Local copy of XXX
%   chanfield           The field whose channels we want to be able to select
%
% RETURNS:
%   hs                  A structure with one fieldname:
%                          'selected_chan_popup', 
%                       which contains a handle to the UI control popup

    global NARFGUI;

    pos = get(parent_handle, 'Position');
    w = pos(3) - 10;
    h = pos(4) - 10;
    hs = [];
    
    m = stack{end}{1};
    mod_idx = length(stack);
    x = xxx{end};

    % If field wasn't defined, use 'output' as the default
    if ~exist('chanfield', 'var'),
        chanfield = m.output;
    end
    
    % Create a channel selector
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'Chan:', ...
        'Units', 'pixels', 'Position', [5 (h-25) 50 25]);
    hs.selected_chan_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [45 (h-25) w-50 25], ...
        'Callback', @(a,b,c) selected_chan_popup_callback());

    function update_channel_popup(stack)
        
        [sf, ~, ~] = get_baphy_plot_controls(stack);
        
        if isfield(x.dat, sf)
            [~, ~, C] = size(x.dat.(sf).(chanfield));
            d = cell(1,C);
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
        hgfeval(get(NARFGUI{mod_idx}.plot_popup, 'Callback'), mod_idx, []);
        drawnow;
    end

    update_channel_popup(stack);

end
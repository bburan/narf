function hs = create_filter_selector_gui(parent_handle, stack, xxx, n_filts)
% hs = create_filter_selector_gui(parent_handle, stack, xxx, n_filts)
%
% Creates a UICONTROL popup menu in parent_handle under the fieldname
% 'selected_filter_popup'. This can then be used to find which filter you
% wish to display for gammatone, elliptic, or other filterbanks.
%
% ARGUMENTS:
%   parent_handle       The panel handle in which the popup will be placed
%
%   stack               Local copy of STACK
%   xxx                 Local copy of XXX
%   n_filts             The number of filters you want to be able to select
%                       between
%
% RETURNS:
%   hs                  A structure with one fieldname:
%                          'selected_filter_popup', 
%                       which contains a handle to the UI control popup
%
    pos = get(parent_handle, 'Position');
    w = pos(3) - 10;
    h = pos(4) - 10;
    hs = [];

    mdl = stack{end};
    mod_idx = length(stack);
    x = xxx{end};
    
    % Create a popup which selects
    uicontrol('Parent', parent_handle, 'Style', 'text', 'Enable', 'on', ...
        'HorizontalAlignment', 'left',  'String', 'Filter#:', ...
        'Units', 'pixels', 'Position', [5 (h-25) w 25]);
    hs.selected_filter_popup = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', 'String', 'NONE', ...
        'Units', 'pixels', 'Position', [5 (h-50) w 25], ...
        'Callback', @(a,b,c) selected_filter_popup_callback());
    
    % Fill that popup with the number of filters
    d = cell(1, n_filts);
    for ii = 1:n_filts
        d{ii} = sprintf('%d',ii);
    end
    set(hs.selected_filter_popup, 'String', char(d));
    set(hs.selected_filter_popup, 'Value', 1);
    
    function selected_filter_popup_callback()
        % Call the plot function again via a sneaky, undocumented callback
        hgfeval(get(mdl.gh.plot_popup,'Callback'), mod_idx, []);
        drawnow;
    end
end
function handles = create_fn_panel(x, y, w, h, parent_handle, state)
% Creates a function panel, with a dropdown to select the function, a data
% table, a plot type dropdown, and a plot axes.
%
% INPUTS:
%    w is the width of the panel created
%    h is the height of the panel created
%    x is the horizontal location is the parent panel, where 0 is the left
%    y is the vertical location, where 0 is the bottom
%    parent_handle is the parent panel this will be created inside
%    state is a global structure which reflects GUI state. IT MUST:
%       1. Have the following fields:
%               pretty_name   A user-readable string
%               plot_name     A user-readable string
%               button_name   A user-readable string
%
% RETURNS: 
%    handles     A structure containing all relevant handles for this

% Create a panel that will hold all the controls for the fn.
handles.fn_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [x y 495 h]);

% Create the popup inside that panel
handles.fn_popup = uicontrol('Parent', handles.fn_panel, ...
    'Style', 'popupmenu', 'Enable', 'on', ...
    'String', state.pretty_name, ...
    'Units', 'pixels', 'Position', [5 (h-30) 300 25], ...
    'Callback', @(a,b,c) disp('FN popup callback TODO.'));

% Create the data table inside that panel
handles.fn_data_table = uitable('Parent', handles.fn_panel, ...
    'Enable', 'on', ...
    'Units', 'pixels', 'Position', [5 35 300 (h-70)], ...
    'CellEditCallback', @(a,b,c) disp('FN Data Table callback TODO.'));

% Create a button at the bottom to apply the function
handles.fn_apply_button = uicontrol('Parent', handles.fn_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', state.button_name, ...
    'Units', 'pixels', 'Position', [185 5 120 25], ...
    'Callback', @(a,b,c) disp('FN Data Table callback TODO.'));

% Create the popup to select which graph to display
handles.plot_type_popup = uicontrol('Parent', handles.fn_panel, ...
    'Style', 'popupmenu', 'Enable', 'on', ...
    'String', state.plot_name, ...
    'Units', 'pixels', 'Position', [310 (h-30) 180 25], ...
    'Callback', @(a,b,c) disp('FN popup callback TODO.'));

% Create a plot-type modification panel for changing visualizations
handles.plot_settings_pane = uipanel('Parent', handles.fn_panel, ...
    'Units', 'pixels', 'Position', [310 5 180 (h-40)]);

% Create the plot axes
% NOTE: Unfortunately, I cannot figure out how to propogate changes from
% the pane to the axes object, which means that if handles.fn_panel is the
% parent of this axes object, then it is not refreshed when parent_handle
% is scrolled. I guess nested update calls isn't working in MATLAB for some
% reason? The ugly workaround is to use parent_handle directly, and
% remember that when you move the fn_panel around, you also need to update
% the axes object associated with that panel.
handles.plot_axes = axes('Parent', parent_handle, ...
            'Units','pixels', 'Position', [520+x 20+y (w-520) (h-25)]);

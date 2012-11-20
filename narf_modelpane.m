function fn_blocks = narf_modelpane(parent_handle, STACK, X)
% A dynamic GUI for chaining together a large number of function calls with
% associated plotting functions to display the input/outputs of the
% function call. The number of functions that may be chained is dynamic and
% can be adjusted on the fly.

fn_blocks = {}; % All handles will go in here

pos = get(parent_handle, 'Position');
w = pos(3);
h = pos(4);

bh = 35;   % Add/del function block button height, +5 pixel padding per side
ph = 170;  % Function panel height

% Create a panel inside it which will slide
handles.container_panel = uipanel('Parent',parent_handle, ...
    'Units','pixels', 'Position',[0 0 w-20 h]);

% Create the scroll bar
handles.container_slider = uicontrol('Parent',parent_handle, ...
    'Style','slider', 'Enable','off', ...
    'Units','pixels', 'Position',[w-20 0 20 h], ...
    'Min',0-eps, 'Max',0, 'Value',0, ...
    'Callback',{@on_scrollbar_slide, handles.container_panel});

function on_scrollbar_slide(hSld, ev, hPan)
    offset = get(hSld,'Value');
    p = get(hPan, 'Position');
    set(hPan, 'Position',[p(1) -offset p(3) p(4)]);
end

function scroll_view_to_bottom(N)
    % Update the size of the parent panel and slider, shift to show bottom
    set(handles.container_panel, 'Position', [0 0 w (bh+ph*N)])
    if (bh+ph*N-h) > 0
        set(handles.container_slider, 'Min', 0, 'Enable', 'on', ...
            'Max', bh+ph*N-h, 'Value', 0); 
    end
    % Use an evil, undocumented function to trigger a callback
    hgfeval(get(handles.container_slider,'Callback'), handles.container_slider, []);
    drawnow;
end

% Define a function that makes new function blocks
function block_handles = create_fn_block_panel(x, y, w, h, parent_handle, state)
    % Creates a function panel, with a dropdown to select the function, a data
    % table, a plot type dropdown, and a plot axes.
    %
    % INPUTS:
    %    pos is the position of the panel [x y w h] where
    %       x is the horizontal location is the parent panel and 0=left
    %       y is the vertical location and 0=bottom
    %       w is the width of the panel created
    %       h is the height of the panel created
    
    %    parent_handle is the parent panel this will be created inside
    %    state is a global structure which reflects GUI state. IT MUST:
    %       1. Have the following fields:
    %               pretty_name   A user-readable string
    %               plot_name     A user-readable string
    %               button_name   A user-readable string
    %
    % RETURNS:
    %    block_handles     A structure containing all relevant block_handles for this
    
    % Create a panel that will hold all the controls for the fn.
    block_handles.fn_panel = uipanel('Parent', parent_handle, ...
        'Units','pixels', 'Position', [x y 495 h]);
    
    % Create the popup inside that panel
    block_handles.fn_popup = uicontrol('Parent', block_handles.fn_panel, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', state.pretty_name, ...
        'Units', 'pixels', 'Position', [5 (h-30) 300 25], ...
        'Callback', @(a,b,c) disp('FN popup callback TODO.'));
    
    % Create the data table inside that panel
    block_handles.fn_data_table = uitable('Parent', block_handles.fn_panel, ...
        'Enable', 'on', ...
        'Units', 'pixels', 'Position', [5 35 300 (h-70)], ...
        'CellEditCallback', @(a,b,c) disp('FN Data Table callback TODO.'));
    
    % Create a button at the bottom to apply the function
    block_handles.fn_apply_button = uicontrol('Parent', block_handles.fn_panel, ...
        'Style', 'pushbutton', 'Enable', 'on', ...
        'String', state.button_name, ...
        'Units', 'pixels', 'Position', [185 5 120 25], ...
        'Callback', @(a,b,c) disp('FN Data Table callback TODO.'));
    
    % Create the popup to select which graph to display
    block_handles.plot_type_popup = uicontrol('Parent', block_handles.fn_panel, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', state.plot_name, ...
        'Units', 'pixels', 'Position', [310 (h-30) 180 25], ...
        'Callback', @(a,b,c) disp('FN popup callback TODO.'));
    
    % Create a plot-type modification panel for changing visualizations
    block_handles.plot_settings_pane = uipanel('Parent', block_handles.fn_panel, ...
        'Units', 'pixels', 'Position', [310 5 180 (h-40)]);
    
    % Create the plot axes
    % NOTE: Unfortunately, I cannot figure out how to propogate changes from
    % the pane to the axes object, which means that if block_handles.fn_panel is the
    % parent of this axes object, then it is not refreshed when parent_handle
    % is scrolled. I guess nested update calls isn't working in MATLAB for some
    % reason? The ugly workaround is to use parent_handle directly, and
    % remember that when you move the fn_panel around, you also need to update
    % the axes object associated with that panel.
    block_handles.plot_axes = axes('Parent', parent_handle, ...
        'Units','pixels', 'Position', [520+x 20+y (w-520) (h-25)]);
end


% Define callbacks for adding and deleting function blocks
function add_fn_block(hObj, evts, unknown)
    n_fns = length(fn_blocks);
    N = n_fns+1;
    
    % Slide all existing panels into their new positions
    for i = 1:n_fns
        set(fn_blocks{i}.plot_axes, 'Position', [520 20+bh+ph*(N-i) (w-520) (ph-25)]);
        set(fn_blocks{i}.fn_panel, 'Position', [0, bh+ph*(N-i) 495 bh]);
    end
    
    % Add a new function block
    menu = [];
    menu.pretty_name = 'Select Function';
    menu.plot_name = 'Select Plot Type';
    menu.button_name = ['Apply Function' num2str(N)];
    fn_blocks{N} = create_fn_block_panel(0, bh, w-20, ph, handles.container_panel, menu);
    
    scroll_view_to_bottom(N);
end

% Define callbacks for adding and deleting function blocks
function del_fn_block(hObj, evts, unknown)
    n_fns = length(fn_blocks);
    
    % Slide all existing panels into their new positions
    for i = 1:n_fns-1
        set(fn_blocks{i}.plot_axes, 'Position', [520 20+bh+ph*(n_fns-i-1) (w-520) (ph-25)]);
        set(fn_blocks{i}.fn_panel,  'Position', [0, bh+ph*(n_fns-i-1) 495 bh]);
    end
    
    % Remove last function block explicitly since we don't trust matlab to
    % garbage collect everything properly
    delete(fn_blocks{n_fns}.plot_axes);
    delete(fn_blocks{n_fns}.plot_settings_pane);
    delete(fn_blocks{n_fns}.fn_apply_button);
    delete(fn_blocks{n_fns}.fn_data_table);
    delete(fn_blocks{n_fns}.fn_popup);
    delete(fn_blocks{n_fns}.fn_panel);
    fn_blocks = fn_blocks(1:n_fns-1);
    scroll_view_to_bottom(n_fns-1);
end

% Create the add/remove function buttons
handles.add_fn_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Add function block', ...
    'Units', 'pixels', 'Position', [5 5 200 25], ...
    'Callback', @add_fn_block);

handles.del_fn_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Delete function block', ...
    'Units', 'pixels', 'Position', [210 5 200 25], ...
    'Callback', @del_fn_block);

% Make the scroll bar dynamically update while being dragged
% hJScrollBar = findjobj(handles.container_slider);
% hJScrollBar.AdjustmentValueChangedCallback = @on_scrollbar_slide;

% TODO: Scroll down once to make sure things are initialized 

% Use an evil, undocumented function to trigger the firstcallback
hgfeval(get(handles.container_slider,'Callback'), handles.container_slider, []);
drawnow;

% COMMENT:
% I would love to make mouse wheel scrolling work for the panel, and not 
% just for the scroll bar, but unfortunately findjobj() cannot return a
% java object for a panel because...MATLAB doesn't use Java swing panels! 
% Therefore, the closest we could do would be to make mouse wheel scrolling
% work just for the scrollbar, although that isn't nearly as much fun.
%
% hJ = findjobj(handles.container_slider);
% hJ.MouseWheelMovedCallback = @(ch, evt, z) disp(get(evt, 'wheelRotation'));

end
function refresh_modelpane = narf_modelpane(parent_handle, modules)
% A dynamic GUI for chaining together a large number of function calls with
% associated plotting functions to display the input/outputs of the
% function call. The number of functions that may be chained is dynamic and
% can be adjusted on the fly. These functions are called 'modules' and have
% an interface that is defined in the documentation. 
%
% ARGUMENTS:
%   parent_handle    The handle to draw the modelpane in
%   modules          A struct of available modules
%
% RETURNS:
%   refresh_fn       A function to refresh the GUI if STACK or XXX have
%                    been edited by another process. 
%
% Once created, the GUI will be able to communicate with other processes by
% modifying global values STACK and XXX.
%
% Any 'gui handles' created here are placed at the right index of the STACK
% data structure under field 'gh'.
% 
% Whenever you modify STACK or XXX without using this GUI, please call
% refresh_modelpane() to update this GUI or weird things may occur.
%
% Unfortunately, limitations of MATLAB seem to prevent the easy way of
% rendering the module blocks, which would be to start at 0 and render each
% module along the negative y axis. If you actually try this, it seems that
% the matlab 'axes' function doesn't render plots placed at negative
% positions. This somewhat obfuscated the implementation of the GUI, since
% everything then has to be shifted around in the positive spaces as the
% number of module blocks changes.

global STACK XXX;

pos = get(parent_handle, 'Position');
w = pos(3); % Width of the parent panel
h = pos(4); % Height of the parent panel

bh = 35;   % Add/del function block button height, +5 pixel padding per side
ph = 170;  % The height of all module panels

% Create a panel inside it which will slide
handles.container_panel = uipanel('Parent',parent_handle, ...
    'Units','pixels', 'Position', [0 0 w-20 h]);

% Create the scroll bar
handles.container_slider = uicontrol('Parent',parent_handle, ...
    'Style','slider', 'Enable','off', ...
    'Units','pixels', 'Position', [w-20 0 20 h], ...
    'Min', 0-eps, 'Max', 0, 'Value', 0, ...
    'Callback', @(h,evts,hds) update_scrolling_on_scrollbar_slide(h, handles.container_panel));

% Create GUIs for any modules that already exist in STACK
% TODO:

% Use as a callback which moves the internal panel around
function update_scrolling_on_scrollbar_slide(hSld, hPan)
    offset = get(hSld,'Value');
    p = get(hPan, 'Position');
    set(hPan, 'Position', [p(1) -offset p(3) p(4)]);
    N = length(STACK);
    top = (N)*ph + bh;
    for yi = 1:N
        set(STACK{yi}.gh.fn_panel, 'Position', [0 (top-ph*yi-offset) 495 bh]);
        set(STACK{yi}.gh.plot_axes, 'Position', [520 (top-ph*yi+20-offset) (w-520) (ph-25)]);
    end
    drawnow;
end

% Make the scroll bar dynamically update while being dragged
hJScrollBar = findjobj(handles.container_slider);
hJScrollBar.AdjustmentValueChangedCallback = @(h, e, v) update_scrolling_on_scrollbar_slide(handles.container_slider, handles.container_panel);

% Update scrollbar and panel sizes whenever the number of modules changes
function update_scrolling_on_module_count_change()
   N = length(STACK);
   hSld = handles.container_slider;
   hPan = handles.container_panel;
   set(hPan, 'Position', [0 0 w (bh+ph*N)]);
   H = max(eps, bh+ph*N-h);
   if H ~= eps
       set(hSld, 'Enable', 'on', 'Min', 0, 'Max', H, 'Value', 0); 
   else
       set(hSld, 'Enable', 'off', 'Min', 0, 'Max', H, 'Value', 0);
   end       
   hgfeval(get(hSld, 'Callback'), hSld, []);
   drawnow;
end

% Called to move the module panels around after a change in STACK's length
function update_module_panel_positions()
    % Calculate the topmost point of the GUI
    N = length(STACK);
    top = (N)*ph + bh;
    for yi = 1:N
        set(STACK{yi}.gh.fn_panel, 'Position', [0 (top-ph*yi) 495 bh]);
        set(STACK{yi}.gh.plot_axes, 'Position', [520 (top-ph*yi+20) (w-520) (ph-25)]);
    end
    drawnow;
end

% Fill a popup menu with modules that are ready to execute at stack depth
% mod_idx
function populate_popup_with_ready_modules(hObject, mod_idx)
    ready_mods = {};
    fns = fieldnames(modules);
    j=1;
    for idx = 1:length(fns);
        if modules.(fns{idx}).isready_pred(STACK(1:mod_idx), XXX(1:mod_idx-1))
            ready_mods{j} = modules.(fns{idx}).pretty_name;
            j = j+1;
        end 
    end
    
    % TODO: Select the first option? Or leave it as previously initialized?
    
    if j == 1
        set(hObject, 'String', 'Select Module');
    else
        set(hObject, 'String', char(ready_mods{:}));
    end
end

% Given a module panel handle and selected module name, update the module
% panel to reflect the available plot options.
function populate_popup_with_available_plots(hObject, mod_name)  
    plot_fns = {};
    fns = fieldnames(modules.(mod_name).plot_fns);
    j = 1;
    for idx = 1:length(fns);
        plot_fns{j} = mod.(fns{idx}).plot_pretty_name;
        j = j+1;
    end
    set(hObject, 'String', char(plot_fns{:}));
end

function module_selected_callback (hObject, mod_idx)
    % Which module was selected?
    c = cellstr(get(hObject,'String'));
    pretty_name = c{get(hObject,'Value')};

    % When my revolution comes, I will destroy languages with for loops and
    % repopulate the earth with higher level functional idioms
    f = [];
    fns = fieldnames(modules);
    for i = 1:length(fns)  
        if isequal(pretty_name, modules.(fns{i}).pretty_name)
            f = modules.(fns{i});
        end    
    end
    
    % If f is not found throw an error
    if ~isequal(f, []) 
        selected = f.fn_name;
    elseif isequal(pretty_name, 'Select Module')
        populate_popup_with_ready_modules(hObject, mod_idx);
        return; 
    else
    	error('Somehow, the selected field name was not found!?');
    end
    
    % All the above work just to find out what was selected! Oof!
    % Set the stack at this point to be the selected module. 
    STACK{mod_idx} = merge_structs(STACK{mod_idx}, modules.(selected));
   
    % Invalidate stack block BELOW this one so errors don't propogate
    % Be sure to destroy the GUI objects too!
    % STACK = STACK(1:mod_idx);
    % TODO:
    
    % Update the data table values
    mhl = modules.(f);
    generic_checkbox_data_table(stack{mod_idx}.gh.fn_table, mdl, mdl.editable_fields); 
    
    % Update the plotting popup menu
    % TODO: populate_popup_with_available_plots();
    
    % Call the default plotting method
    % TODO: 
end

% Define a function that makes new module gui blocks
function gh = create_mod_block_panel(parent_handle, mod_idx)
    % Returns a module panel placed at the origin. 
    %
    % The returned 'gui handle' gh has these fields:
    %   .fn_panel    The panel which encapsulates all active controls.
    %   .plot_axes   The axes object which is used for plotting.
    %   .fn_popup    Popup menu to select the module.
    %   .fn_table    A data table to editing the module parameters.
    %   .plot_popup  Popup menu to select which plot to view.
    %   .plot_panel  A panel for user-defined plot controls.
    %   .fn_apply    A button to apply the function.
    % 
    %
    %    parent_handle is the parent panel this will be created inside
    %    state is a global structure which reflects GUI state. IT MUST:
    %       1. Have the following fields:
    %               pretty_name   A user-readable string
    %               plot_name     A user-readable string
    %               button_name   A user-readable string
    %
    % RETURNS:
    %    gh     A structure containing all relevant handles
    
    % Create a panel that will hold all the controls for the fn.
    gh.fn_panel = uipanel('Parent', parent_handle, ...
        'Units','pixels', 'Position', [0 0 495 ph]);
    
    % Create the popup to select which graph to display
    gh.plot_popup = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'Uninitialized', ...
        'Units', 'pixels', 'Position', [310 (ph-30) 180 25], ...
        'Callback', @(a,b,c) fprintf('plot_popup callback for block %d\n', mod_idx));
    
    % Create a plot-type modification panel for changing visualizations
    gh.plot_panel = uipanel('Parent', gh.fn_panel, ...
        'Units', 'pixels', 'Position', [310 5 180 (ph-40)]);
    
    % Create the popup inside that panel
    gh.fn_popup = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'Select Module', ...
        'Units', 'pixels', 'Position', [5 (ph-30) 300 25], ...
        'Callback', @(h, evts, hdls) module_selected_callback(h, mod_idx));

    % Create the data table inside that panel
    gh.fn_table = uitable('Parent', gh.fn_panel, ...
        'Enable', 'on', ...
        'Units', 'pixels', 'Position', [5 35 300 (ph-70)], ...
        'CellEditCallback', @(a,b,c) disp('FN Data Table callback TODO.\n'));
    
    % Create a button at the bottom to apply the function
    gh.fn_apply = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'pushbutton', 'Enable', 'on', ...
        'String', 'Apply!', ...
        'Units', 'pixels', 'Position', [185 5 120 25], ...
        'Callback', @(a,b,c) fprintf('fn_apply callback for mod %d\n', mod_idx));
    
    % Create the plot axes
    % NOTE: Unfortunately, I cannot figure out how to propogate changes from
    % the pane to the axes object, which means that if gh.fn_panel is the
    % parent of this axes object, then it is not refreshed when parent_handle
    % is scrolled. I guess nested update calls isn't working in MATLAB for some
    % reason? The ugly workaround is to use parent_handle directly, and
    % remember that when you move the fn_panel around, you also need to update
    % the axes object associated with that panel.
    gh.plot_axes = axes('Parent', parent_handle, ...
        'Units','pixels', 'Position', [520 20 (w-520) (ph-25)]);
end

% Callback for creating a module block
function add_mod_block()
    % Add a new model block, update the stack, and finally the view
    idx = (length(STACK)+1);
    STACK{idx}.gh = create_mod_block_panel(parent_handle, idx);
    update_module_panel_positions();    
    update_scrolling_on_module_count_change();
end

% Callback for deleting a module block
function del_mod_block()
    n_fns = length(STACK);
    
    if n_fns == 0
        return;
    end
    
    % Remove GUI handles explicitly since I don't trust matlab to
    % garbage collect everything properly
    delete(STACK{n_fns}.gh.plot_axes);
    delete(STACK{n_fns}.gh.plot_popup);
    delete(STACK{n_fns}.gh.plot_panel);
    delete(STACK{n_fns}.gh.fn_apply);
    delete(STACK{n_fns}.gh.fn_table);
    delete(STACK{n_fns}.gh.fn_popup);
    delete(STACK{n_fns}.gh.fn_panel);
    
    % Update the stack, the remaining module posiitons, and the view
    STACK = STACK(1:n_fns-1);
    update_module_panel_positions();    
    update_scrolling_on_module_count_change();
end

% Create the add/remove function buttons
handles.add_fn_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Add function block', ...
    'Units', 'pixels', 'Position', [5 5 200 25], ...
    'Callback', @(h,b,c) add_mod_block());

handles.del_fn_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Delete function block', ...
    'Units', 'pixels', 'Position', [210 5 200 25], ...
    'Callback', @(h,b,c) del_mod_block());

% Use an evil, undocumented function to trigger the firstcallback
hgfeval(get(handles.container_slider,'Callback'), handles.container_slider, []);
drawnow;

% COMMENT:
% I would love to make mouse wheel scrolling work with the cursor over the
% entire panel, and not when the cursor is over  the scroll bar. 
% Unfortunately findjobj() cannot return a  java object for a panel 
% because...MATLAB doesn't use Java swing panels! 
% Therefore, the closest we could do would be to make mouse wheel scrolling
% work just for the scrollbar, although that isn't nearly as much fun.
%
% hJ = findjobj(handles.container_slider);
% hJ.MouseWheelMovedCallback = @(ch, evt, z) disp(get(evt, 'wheelRotation'));

end
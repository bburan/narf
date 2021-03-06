function refresh_modelpane = narf_modelpane(parent_handle)
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

global STACK XXX META NARFGUI MODULES;

MODULES = scan_directory_for_modules();

if ~exist('parent_handle','var')
    parent_handle = figure('Menubar','figure', 'Resize','off', 'MenuBar', 'none', ...              
             'Units','pixels', 'Position', [20 50 1300 max(800, min(170*length(STACK)+40, 1100))]);
    % 'ResizeFcn', @update_windowsize, ... Someday write this!
end

if isprop(parent_handle, 'Menubar')
    set(parent_handle, 'Resize','off', 'Menubar', 'none', ...
       'Name', 'NARF Modelpane', 'NumberTitle', 'off');
end
   
pos = get(parent_handle, 'Position');
w = pos(3); % Width of the parent panel
h = pos(4); % Height of the parent panel

bh = 65;   % Add/del function block button height, +5 pixel padding per side
ph = 170;  % The height of all module panels

% Create a panel inside it which will slide
handles.container_panel = uipanel('Parent',parent_handle, ...
    'Units','pixels', 'Position', [0 0 w-20 h]);

% Create the scroll bar
handles.container_slider = uicontrol('Parent',parent_handle, ...
    'Style','slider', 'Enable','off', 'SliderStep',[0.2 1.0], ...
    'Units','pixels', 'Position', [w-20 0 20 h], ...
    'Min', 0-eps, 'Max', 0, 'Value', 0, ...
    'Callback', @(h,evts,hds) update_panel_positions());

% Move all panels around to their proper positions
function update_panel_positions()
    hSld = handles.container_slider;
    hPan = handles.container_panel;
    offset = get(hSld,'Value');
    p = get(hPan, 'Position');
    set(hPan, 'Position', [p(1) -offset p(3) p(4)]);
    N = length(NARFGUI);
    top = (N)*ph + bh;
    for yi = 1:N
        p = get(NARFGUI{yi}.fn_panel, 'Position');
        set(NARFGUI{yi}.fn_panel, 'Position', [0 (top-ph*yi-offset) p(3) p(4)]);
        p = get(NARFGUI{yi}.plot_axes, 'Position');
        set(NARFGUI{yi}.plot_axes, 'Position', [550 (top-ph*yi+20-offset) p(3) p(4)]);
    end
    drawnow;
end

% Update scrollbar and panel sizes whenever the number of modules changes
function update_scrollbar_size()
   N = length(NARFGUI);
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

% Fill a popup menu with modules that are ready at stack depth 'mod_idx'
function update_ready_modules(mod_idx)
    ready_mods = {'Select Module'};
    fns = fieldnames(MODULES);
    j=2;
    if mod_idx > length(XXX)  % If there is no data, don't do anything
        return
    end
    for idx = 1:length(fns);
        if MODULES.(fns{idx}).isready_pred(STACK, XXX)
            ready_mods{j} = MODULES.(fns{idx}).pretty_name;
            j = j+1;
        end 
    end
    
    % TODO: Select the first option? Or leave it as previously initialized?
    set(NARFGUI{mod_idx}.fn_popup, 'String', char(ready_mods{:}));
end

% Given a module panel handle and selected module name, update the module
% panel to reflect the available plot options.
function update_available_plots(mod_idx)  
    plot_fns = {};
    pfns = STACK{mod_idx}{1}.plot_fns;
    for idx = 1:length(pfns);
        plot_fns{idx} = pfns{idx}.pretty_name;
    end
    if isempty(plot_fns)
        set(NARFGUI{mod_idx}.plot_popup, 'String', 'Unavailable'); 
    else
        set(NARFGUI{mod_idx}.plot_popup, 'String', char(plot_fns{:}));
    end
end

% A method for updating a simple fieldname/value uitable
function update_uitable(mytable, mystruct, myfields)
    l = length(myfields);
    c = cell(l,2);
    for i = 1:l
        if ~isfield(mystruct, myfields{i})
            log_err('Could not find field: %s', myfields{i});
        end
        c{i,1} = myfields{i};
        c{i,2} = write_readably(mystruct.(myfields{i})); % Ensure data becomes a str
    end
    set(mytable, 'Data', c);
    drawnow;
end

% A method for abstractly updating a data table in a module
function update_checkbox_uitable(mytable, mystruct, myfields)
    l = length(myfields);
    c = cell(l,3);
    for i = 1:l
        if ~isfield(mystruct, myfields{i})
            error('Could not find field: %s', myfields{i});
        end
        if isfield(mystruct, 'fit_fields')
            c{i,1} = any(strcmp(myfields{i}, mystruct.fit_fields));
        else
            c{i,1} = false;
        end
        c{i,2} = myfields{i};
        c{i,3} = write_readably(mystruct.(myfields{i})); % Ensure data becomes a str
    end
    set(mytable, 'ColumnName', {'Fit?', 'Field', 'Value'});
    set(mytable, 'ColumnEditable', [true false true]);
    set(mytable, 'ColumnWidth', {25 150 100});
    set(mytable, 'RowName', {});
    set(mytable, 'Data', c);
    
    drawnow;
end

function module_selected_callback(mod_idx)
    % Which module was selected?
    c = cellstr(get(NARFGUI{mod_idx}.fn_popup,'String'));
    pretty_name = c{get(NARFGUI{mod_idx}.fn_popup,'Value')};
    
    % When my revolution comes, I will destroy languages with for loops and
    % repopulate the earth with higher level functional idioms
    f = [];
    fns = fieldnames(MODULES);
    for i = 1:length(fns)  
        if isequal(pretty_name, MODULES.(fns{i}).pretty_name)
            f = MODULES.(fns{i});
        end    
    end
    
    % If f is not found throw an error
    if ~isequal(f, []) 
        selected = f.name;
    elseif isequal(pretty_name, 'Select Module')
        update_ready_modules(mod_idx);
        return; 
    else
    	error('Somehow, the selected field name was not found!?');
    end
    
    % All the above work just to find out what was selected! Oof!
    % Set the stack at this point to be the selected module. 
    STACK{mod_idx} = {MODULES.(selected)};

    % Invalidate data beyond this point
    XXX = XXX(1:mod_idx);
    
    % Update the data table values 
    m = STACK{mod_idx}{1};
    update_checkbox_uitable(NARFGUI{mod_idx}.fn_table, m, m.editable_fields); 
    
    % Update the plotting popup menu, but leave it disabled
    update_available_plots(mod_idx);
    set(NARFGUI{mod_idx}.plot_popup, 'Enable', 'off');   
    set(NARFGUI{mod_idx}.fn_recalc, 'Value', false);   
    set(NARFGUI{mod_idx}.fn_replot, 'Value', false);   
        
end


% When the user presses the button, apply the function
function module_apply_callback(mod_idx)
    m = STACK{mod_idx}{1};
    XXX = XXX(1:mod_idx);  % Invalidate later data so it cannot be 
                           % accidentally used by this or later functions
    
    % Check again that we are able to run; the user may have changed an
    % upstream part of the stack, leaving this function actually invalid
    if m.isready_pred(STACK(1:mod_idx), XXX)
        % Apply the function
        calc_xxx(mod_idx, mod_idx);
        % Enable graphing
        set(NARFGUI{mod_idx}.plot_popup, 'Enable', 'on');
        % Build the plot panel, now that we know XXX{mod_idx+1}
        % Only build if it doesn't already exist
        if isfield(m, 'plot_gui_create_fn') && ~isfield(NARFGUI{mod_idx}, 'plot_gui')
            NARFGUI{mod_idx}.plot_gui = m.plot_gui_create_fn(NARFGUI{mod_idx}.plot_panel, STACK(1:mod_idx), XXX(1:mod_idx+1));
        end
        % Trigger a redraw of the gui
        hgfeval(get(NARFGUI{mod_idx}.plot_popup,'Callback'), mod_idx, []);
        drawnow;

         % If auto-calc of the NEXT fn is checked, go down the stack
        if length(STACK) > mod_idx + 1 && ...
           isequal(true, get(NARFGUI{mod_idx+1}.fn_recalc, 'Value'))
           module_apply_callback(mod_idx+1);
        end
    else
        fprintf('Sorry, Dave, I''m afraid I can''t do that...\n');
    end
end

% When the user changes the popup selection update the plot
% TODO: Right now this is used as the general-use interface to the
% plotting subsystem and it's a little hacky. 
function module_plot_callback(mod_idx)
    m = STACK{mod_idx}{1};

    % Get the selected plot function
    idx = get(NARFGUI{mod_idx}.plot_popup, 'Value');    
    
    % If the the index is a valid one and the function has been run (which
    % would make XXX(mod_idx+1) useful to us), we can plot.
    if (idx > 0 && idx <= length(m.plot_fns) && length(XXX) >= mod_idx+1)
        % Set the axes, clear it, and run the plot function
        axes(NARFGUI{mod_idx}.plot_axes);
        cla; legend off;
        sel = narfgui_widget_selected_values(NARFGUI(1:mod_idx));        
        m.plot_fns{idx}.fn(sel, STACK(1:mod_idx), XXX(1:mod_idx+1));
        replot_from_depth(mod_idx+1);
    end 
    
    % Special case: If a module has no plotting, pass it along.
    if isempty(m.plot_fns)
        replot_from_depth(mod_idx+1);
    end
    
end

function recalc_from_depth(mod_idx)
    % Like calc_xxx() but for use as a gui callback
    % Try to rebuild XXX, starting at stack depth mod_idx
    % Unlike calc_xxx(), it stops when a module's autoapply isn't checked
    for ii = mod_idx:length(STACK);
        if isfield(NARFGUI{ii}, 'fn_recalc') && get(NARFGUI{ii}.fn_recalc, 'Value')
            module_apply_callback(ii);
        else
            return
        end
    end
end

function replot_from_depth(mod_idx)
    % Try to replot everything from mod_idx onward
    % Stop trying to plot if you reach an unchecked 'autoplot' checkbox
    if mod_idx > length(STACK)
        return
    end
    if isfield(NARFGUI{mod_idx}, 'fn_replot') && get(NARFGUI{mod_idx}.fn_replot, 'Value')
        module_plot_callback(mod_idx); % Recurses
    end
end

% When the data table changes, invalidate the data, plot and plot gui
function module_data_table_callback(mod_idx)
    
    % Extract the structure of the field
    s = extract_field_val_pairs(NARFGUI{mod_idx}.fn_table, 2, 3);
    
    if length(STACK{mod_idx}) > 1 &&...
        ~strcmp('Yes', questdlg('Overwrite ALL parameter sets? This GUI has no way of modifying a single parameter set yet.'))
        s = [];
    end
    
    % If there was an error extracting the fields, reset the GUI
    if isempty(s)
        update_checkbox_uitable(NARFGUI{mod_idx}.fn_table, ...
                                STACK{mod_idx}{1}, ...
                                STACK{mod_idx}{1}.editable_fields);
        return;
    end
    
    % Throw an error if you are trying to edit a field that has different
    % values across each parameter set
    % OR, make sure that you are editing only the SELECTED parameter set
    
    for ii = 1:length(STACK{mod_idx})    
        % Insert the field values into the stack
        for fs = fieldnames(s)', fs=fs{1};
            STACK{mod_idx}{ii}.(fs) = s.(fs);
        end
    
        % Update the selected fields
        STACK{mod_idx}{ii}.fit_fields = extract_checked_fields(NARFGUI{mod_idx}.fn_table, 1, 2);
    
        % Give the module a chance to run its own code
        STACK{mod_idx}{ii} = STACK{mod_idx}{ii}.mdl(STACK{mod_idx}{ii});  
    end
    
    % Request a recalculation from this point onwards. 
    recalc_from_depth(mod_idx);
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
    %   .fn_recalc   A checkbox to determine if it should auto-recalc
    %   .fn_replot   A checkbox to determine if it should auto-plot
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
        'Callback', @(a,b,c) module_plot_callback(mod_idx));
    
    % Create a plot-type modification panel for changing visualizations
    gh.plot_panel = uipanel('Parent', gh.fn_panel, ...
        'Units', 'pixels', 'Position', [310 5 180 (ph-40)]);
    
    % Create the popup inside that panel
    gh.fn_popup = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', 'Select Module', ...
        'Units', 'pixels', 'Position', [5 (ph-30) 300 25], ...
        'Callback', @(h, evts, hdls) module_selected_callback(mod_idx));

    % Create the data table inside that panel
    gh.fn_table = uitable('Parent', gh.fn_panel, ...
        'Enable', 'on', ...
        'Units', 'pixels', 'Position', [5 35 300 (ph-70)], ...
        'CellEditCallback', @(a,b,c) module_data_table_callback(mod_idx));
       
    % Create a button at the bottom to apply the function
    gh.fn_apply = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'pushbutton', 'Enable', 'on', ...
        'String', 'Apply!', ...
        'Units', 'pixels', 'Position', [185 5 120 25], ...
        'Callback', @(a,b,c) module_apply_callback(mod_idx));
    
    % Create a automatic recalc checkbox
    gh.fn_recalc = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'checkbox', 'Enable', 'on', 'Value', true, ...
        'String', 'AutoApply', ...
        'Units', 'pixels', 'Position', [5 5 90 25], ...
        'Callback', @(a,b,c) enable_or_disable_button());
    
    function enable_or_disable_button()
        autoapply = get(gh.fn_recalc, 'Value');
        if autoapply
            set(gh.fn_apply, 'Enable', 'off');
            module_apply_callback(mod_idx);
        else
            set(gh.fn_apply, 'Enable', 'on');
        end
    end
    
    % Create an automatic replot checkbox
    gh.fn_replot = uicontrol('Parent', gh.fn_panel, ...
        'Style', 'checkbox', 'Enable', 'on', 'Value', true, ...
        'String', 'AutoPlot', ...
        'Units', 'pixels', 'Position', [95 5 90 25]);
    
    % Create the plot axes
    % NOTE: Unfortunately, I cannot figure out how to propogate changes from
    % the pane to the axes object, which means that if gh.fn_panel is the
    % parent of this axes object, then it is not refreshed when parent_handle
    % is scrolled. I guess nested update calls isn't working in MATLAB for some
    % reason? The ugly workaround is to use parent_handle directly, and
    % remember that when you move the fn_panel around, you also need to update
    % the axes object associated with that panel.
    gh.plot_axes = axes('Parent', parent_handle, ...
        'Units','pixels', 'Position', [550 20 (w-550-25) (ph-25)]);
end

% Callback for creating a module block
function add_mod_block()
    if length(NARFGUI) > length(STACK)
        return
    end
    % Add a new model block, update the stack, and finally the view
    idx = (length(NARFGUI)+1);   
    NARFGUI{idx} = create_mod_block_panel(parent_handle, idx);
    % Trigger the popup callback to initialize it
    hgfeval(get(NARFGUI{idx}.fn_popup, 'Callback'), idx, []);
    % Update the view
    update_scrollbar_size();
    update_panel_positions();      
end

% Callback for deleting a module block
function del_mod_block()
    idx = length(NARFGUI);
    
    if idx == 0 
        return;
    end
    
    delete_module_gui(idx);
    NARFGUI = NARFGUI(1:idx-1);
    
    % Update the stack, the remaining module positions, and the view
    if length(STACK) == idx,
        STACK = STACK(1:idx-1);
        XXX = XXX(1:idx);
    end
    update_scrollbar_size();
    update_panel_positions();    
    
end

function save_model_callback (a, b, c)
    [filename, pathname] = uiputfile({'*.mat'}, ...
                                     'Save Model Stack As', ...
                                     META.modelpath);
    if isequal(filename,0) || isequal(pathname,0) 
        % Canceled
    else
        save_model([pathname filesep filename], STACK, XXX, META);
        db_insert_model();
    end
end

function load_model_callback (a, b, c)
    if isfield(META, 'modelpath')
        [filename, pathname] = uigetfile('*.mat', ...
                                     'Select Model Stack', META.modelpath);    
    else
        [filename, pathname] = uigetfile('*.mat', ...
                                     'Select Model Stack');
    end
    
	if isequal(filename,0) || isequal(pathname,0) 
        % Canceled
    else
        delete_all_module_guis();
        load_model([pathname filesep filename]);
        calc_xxx(1);
        rebuild_gui_from_stack();
    end
end

% Create the add/remove function buttons
handles.add_fn_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Add Module', ...
    'Units', 'pixels', 'Position', [5 5 140 25], ...
    'Callback', @(h,b,c) add_mod_block());

handles.del_fn_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Delete Module', ...
    'Units', 'pixels', 'Position', [5 30 140 25], ...
    'Callback', @(h,b,c) del_mod_block());

% Create save/load model buttons
handles.save_model_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Save Module Stack', ...
    'Units', 'pixels', 'Position', [w-165 5 140 25], ...
    'Callback', @save_model_callback);

handles.load_model_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Load Module Stack', ...
    'Units', 'pixels', 'Position', [w-165 30 140 25], ...
    'Callback', @load_model_callback);

function update_any_changed_tables_and_recalc()
    first_fit_depth = 1;
    % Loop through and update any data tables which changed
    for ii = 1:length(STACK)
        if isfield(STACK{ii}{1}, 'fit_fields')
            first_fit_depth = min(first_fit_depth, ii);
            update_checkbox_uitable(NARFGUI{ii}.fn_table, ...
                STACK{ii}{1}, ...
                STACK{ii}{1}.editable_fields);
        end
    end
    % Recalc from the first fittable point all the way to the end
    calc_xxx(first_fit_depth);
    module_plot_callback(first_fit_depth);
    update_modelinfo_text();
end

handles.fitter_dropdown = uicontrol('Parent', handles.container_panel, ...
    'Style', 'popupmenu', 'Enable', 'on', ...
    'String', 'Uninitialized', ...
    'Units', 'pixels', 'Position', [160 30 150 25]);

handles.fit_button = uicontrol('Parent', handles.container_panel, ...
    'Style', 'pushbutton', 'Enable', 'on', ...
    'String', 'Fit!', ...
    'Units', 'pixels', 'Position', [160 5 150 25], ...
    'Callback', @(a, b, c) wrapper_for_fit());

function fitter_dropdown_init()
	global NARF_FITTERS_PATH;
    fitters = dir2cell([NARF_FITTERS_PATH filesep '*.m']);
    fitters = cellfun(@(x) x(1:end-2), fitters, 'UniformOutput', false);
	set(handles.fitter_dropdown, 'String', char(fitters{:}));
end

fitter_dropdown_init();

function wrapper_for_fit()
    fitter = popup2str(handles.fitter_dropdown);
    func = str2func(fitter);
    func();
    update_any_changed_tables_and_recalc();
end

handles.modelinfo_text = uicontrol('Parent', handles.container_panel, ...
    'Style', 'text', 'Enable', 'on', 'HorizontalAlignment', 'left', ...
    'String', 'Uninitialized', ...
    'Units', 'pixels', 'Position', [320 5 500 bh-15]);

function update_modelinfo_text()
    if isfield(META, 'modelname') 
        set(handles.modelinfo_text, 'String', ...
            sprintf('Batch:     %d\nCellid:    %s\nModel:    %s', ...
                     META.batch, XXX{1}.cellid,  META.modelname));
    end
end

function rebuild_gui_from_stack()
    delete_all_module_guis();
        
    % Create GUIs for any modules that already exist in STACK
    for ii = 1:length(STACK)
        NARFGUI{ii} = create_mod_block_panel(parent_handle, ii);
        % Set the popup to display the selected model type
        % TODO: Also make it display 'Select Module' so it's not locked in
        set(NARFGUI{ii}.fn_popup, 'String', STACK{ii}{1}.pretty_name, 'Value', 1);
        update_checkbox_uitable(NARFGUI{ii}.fn_table, STACK{ii}{1}, ...
                                    STACK{ii}{1}.editable_fields); 
        update_available_plots(ii);                         
    end
    
    % Trigger the first drawing
    update_scrollbar_size();
    update_panel_positions();

    % The plotting and plot_gui's need to know about the XXX struct to work
    % So we rebuild the XXX struct as needed, and then call replot. 
    for ii = 1:length(STACK) 
        m = STACK{ii}{1};
        
        % Recompute the data 
        if length(XXX) <= ii
            calc_xxx(ii, ii+1);
        else
            calc_xxx(ii);
        end
           
        % Delete any existing plot guis        
        if isfield(NARFGUI{ii}, 'plot_gui')
             rmfield(NARFGUI{ii}, 'plot_gui');
        end
    
        % Now everything is ready to rebuild the GUI if there is one
        if isfield(STACK{ii}{1}, 'plot_gui_create_fn')
            NARFGUI{ii}.plot_gui = STACK{ii}{1}.plot_gui_create_fn(NARFGUI{ii}.plot_panel, STACK(1:ii), XXX(1:ii+1));
        end
    end
    
    % If any modules are defined, we should now be able to update the graphs
    % I would have done this earlier, but module_plot_callback works down the
    % whole stack and is recursive, so it can't be run until everything is
    % ready.
    if 0 < length(STACK)
        module_plot_callback(1);
        update_modelinfo_text();
    end   
    
end

rebuild_gui_from_stack();

% Make the scroll bar dynamically update whenever it is being dragged
hJScrollBar = findjobj(handles.container_slider);
hJScrollBar.AdjustmentValueChangedCallback = @(h, e, v) update_panel_positions();

% Define a close window callback which will remove all GUI hooks from the
% STACK. Useful if you close the GUI and want to open it up again without
% recomputing the whole damn STACK.
function delete_gui_and_close(a,b,c)
%     selection = questdlg('Close the GUI? (STACK and XXX will still be there)',...
%        'Close Request Function', 'Yes', 'No', 'Yes'); 
    selection = 'Yes'; 
    switch selection, 
       case 'Yes',
          delete_all_module_guis();
          delete(parent_handle);
       case 'No'
       return 
    end
end

% % Unfortunately, I don't know how to check if parent_handle is a panel or a
% % figure, so I'm using a try/catch instead of a proper solution.
% % if ~isuipanel(parent_handle)
% try 
%     set(parent_handle, 'CloseRequestFcn', @delete_gui_and_close);
% catch
    % Do nothing if it failed
%end

refresh_modelpane = @rebuild_gui_from_stack;

end

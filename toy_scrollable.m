function fn_blocks = toy_scrollable()
% A dynamic GUI for chaining together a large number of function calls with
% associated plotting functions to display the input/outputs of the
% function call. The number of functions that may be chained is dynamic and
% can be adjusted on the fly.

fn_blocks = {}; % All handles

w = 1600;  % Total pane width, including slider
h = 800;   % Total pane height
bh = 35;   % Add/del function block button height, +5 pixel padding per side
ph = 170;  % Function panel height

% Create the figure
handles.hFig = figure('Menubar','figure', 'Resize','off', ...
    'Units','pixels', 'Position',[20 50 w h]);

% Create a panel inside it which will slide
handles.container_panel = uipanel('Parent',handles.hFig, ...
    'Units','pixels', 'Position',[0 0 w-20 h]);

% Create the scroll bar
handles.container_slider = uicontrol('Parent',handles.hFig, ...
    'Style','slider', 'Enable','off', ...
    'Units','pixels', 'Position',[w-20 0 20 h], ...
    'Min',0-eps, 'Max',0, 'Value',0, ...
    'Callback',{@on_scrollbar_slide, handles.container_panel});

% Make the scroll bar dynamically update while being dragged
hJScrollBar = findjobj(handles.container_slider);
hJScrollBar.AdjustmentValueChangedCallback = @(a,b,c) on_scrollbar_slide(handles.container_slider, [], handles.container_panel);

function on_scrollbar_slide(hSld, ev, hPan)
    offset = get(hSld,'Value');
    p = get(hPan, 'Position');
    set(hPan, 'Position',[p(1) -offset p(3) p(4)])
end

function scroll_view_to_bottom(N)
    % Update the size of the parent panel and slider, shift to show bottom
    set(handles.container_panel, 'Position', [0 0 w max(h, (bh+ph*N))])
    set(handles.container_slider, 'Min', 0, 'Enable', 'on', ...
        'Max', max(h, (bh+ph*N)), 'Value', 0); 
    
    % Use an evil, undocumented function to trigger a callback
    hgfeval(get(handles.container_slider,'Callback'), handles.container_slider, []);
    drawnow;
end

% Define a function that makes new function blocks





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
    fn_blocks{N} = create_fn_panel(0, bh, w-20, ph, handles.container_panel, menu);
    
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

% Move later to its own function
% for i = 1:6
%     
%     % Slide the position of all existing panels up by 1. 
%     
%     % Slide the slider bar down to the bottom
%     
%     menu = [];
%     menu.pretty_name = 'Select Function';
%     menu.plot_name = 'Select Plot Type';
%     menu.button_name = ['Apply Function' num2str(i)];
%     fn_blocks{i} = create_fn_panel(0, ph*(i-1), w-20, ph, handles.container_panel, menu);
%     
%     % Update the size of the parent panel and slider, shift to show top
%set(handles.container_panel, 'Position', [0 (ph*i-h) w (ph*i)])
%set(handles.container_slider, 'Min', 0, 'Enable', 'on', ...
        %'Max', (ph*N-h), 'Value', 0); 
%     
% end

% Use an evil, undocumented function to trigger a callback
%hgfeval(get(handles.container_slider,'Callback'), handles.container_slider, []);
%drawnow;

% Make scroll bar work dynamically.

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
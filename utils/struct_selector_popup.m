function options = struct_selector_popup(input_options)
% selected = struct_selector_popup(options)
%
% Given struct THESTRUCT (and optionally a second struct SELECTED), it
% creates a popup that allows you to select and order elements to be
% displayed, which will then be returned as a cell array.

options = input_options;

quitting = false;

lw = 300;
w = 500;
bh= 25;
h = 500;
s = 5;

parent_handle = figure('Menubar','none', 'Resize','off', ...
     'Units','pixels', 'Position', [20 50 w h],...
     'Name', 'Column Selection Pop-up', 'NumberTitle', 'off');
 
tab = uicontrol('Parent', parent_handle, 'style', 'listbox', ...
    'Enable', 'on', 'Units', 'pixels', ...
    'Position', [s s lw-2*s h-2*s]);

button_up = uicontrol('Parent', parent_handle, 'style', 'pushbutton', ...
    'Enable', 'on', 'Units', 'pixels', ...
    'String', 'Move Up', ...
    'Position', [lw+s h-bh*1-s w-lw-2*s bh], ...
    'Callback', @up_callback);

toggle = uicontrol('Parent', parent_handle, 'style', 'pushbutton', ...
    'Enable', 'on', 'Units', 'pixels', ...
    'String', '<Toggle Visibility>', ...
    'Position', [lw+s h-bh*2-s w-lw-2*s bh], ...
    'Callback', @toggle_callback);

button_down = uicontrol('Parent', parent_handle, 'style', 'pushbutton', ...
    'Enable', 'on', 'Units', 'pixels', ...
    'String', 'Move Down', ...
    'Position', [lw+s h-bh*3-s w-lw-2*s bh], ...
    'Callback', @down_callback);

button_reset = uicontrol('Parent', parent_handle, 'style', 'pushbutton', ...
    'Enable', 'on', 'Units', 'pixels', ...
    'String', 'Reset to Defaults', ...
    'Position', [lw+s s+bh w-lw-2*s bh], ...
    'Callback', @reset_callback);

button_done = uicontrol('Parent', parent_handle, 'style', 'pushbutton', ...
    'Enable', 'on', 'Units', 'pixels', ...
    'String', 'Save & Exit', ...
    'Position', [lw+s s w-lw-2*s bh], ...
    'Callback', @done_callback);

function refresh_display()
    d = {};
    for ii = 1:length(options)
        if options{ii}{2}
            d{end+1} = options{ii}{1};
        else     
            d{end+1} = [options{ii}{1} ' (Hidden)'];
        end
    end            
    set(tab, 'String', d);
end

function up_callback(~,~,~)
    idx = get(tab, 'Value');
    if (idx > 1)
        tmp = options{idx-1};
        options{idx-1} = options{idx};
        options{idx} = tmp;
        set(tab, 'Value', idx-1);
    end
    refresh_display();
end

function toggle_callback(~,~,~)
    idx = get(tab, 'Value');
    options{idx}{2} = ~options{idx}{2};   
    refresh_display();
end

function down_callback(~,~,~)
    idx = get(tab, 'Value');
    if (idx < length(options))
        tmp = options{idx+1};
        options{idx+1} = options{idx};
        options{idx} = tmp;       
        set(tab, 'Value', idx+1);
    end
    refresh_display();
end

function reset_callback(~,~,~)
    options = input_options;
    refresh_display();
end

function done_callback(~,~,~)
    quitting = true;
    close(parent_handle);
end

refresh_display();

while ~(quitting) && ishandle(parent_handle)
    pause(0.2);
end

end
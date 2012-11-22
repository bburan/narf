function x = make_x_modifier(x)

    function increment_it()
        x = x + 1;
        fprintf('incremented to: %d\n', x);
    end
   
    pf = figure('Menubar','figure', 'Resize','off', ...
        'Units','pixels', 'Position', [20 50 500 500]);
    
    mybutton = uicontrol('Parent', pf, ...
        'Style', 'pushbutton', 'Enable', 'on', ...
        'String', 'Push to increment', ...
        'Units', 'pixels', 'Position', [185 5 120 25], ...
        'Callback', @(a,b,c) increment_it());
    
end
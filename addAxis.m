function hAx = addAxis(handles)
    %# look for previous axes
    ax = findobj(handles.hPan, 'type','axes');

    if isempty(ax)
        %# create first axis
        hAx = axes('Parent',handles.hPan, ...
            'Units','normalized', 'Position',[0.13 0.11 0.775 0.815]);
        set(hAx, 'Units','pixels');

    else
        %# get height of figure
        p = get(handles.hFig, 'Position');
        h = p(4);

        %# increase panel height, and shift it to show new space
        p = get(handles.hPan, 'Position');
        set(handles.hPan, 'Position',[p(1) p(2)-h p(3) p(4)+h])

        %# compute position of new axis: append on top (y-shifted)
        p = get(ax, 'Position');
        if iscell(p), p = cell2mat(p); end
        p = [p(1,1) max(p(:,2))+h p(1,3) p(1,4)];

        %# create the new axis
        hAx = axes('Parent',handles.hPan, ...
            'Units','pixels', 'Position',p);

        %# adjust slider, and call its callback function
        mx = get(handles.hSld, 'Max');
        set(handles.hSld, 'Max',mx+h, 'Min',0, 'Enable','on')
        %#set(handles.hSld, 'Value',mx+h)       %# scroll to new space
        hgfeval(get(handles.hSld,'Callback'), handles.hSld, []);
    end

    %# force GUI update
    drawnow
end
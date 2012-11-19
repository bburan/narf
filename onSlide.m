function onSlide(hSld, ev, hPan)
    %# slider value    
    offset = get(hSld,'Value');

    %# update panel position
    p = get(hPan, 'Position');  %# panel current position
    set(hPan, 'Position',[p(1) -offset p(3) p(4)])
    
    % Run callbacks on the axes, since they are not JAVA and not being
    % updated properly automatically.     
%     for i = 1:length(my_handles)
%         hgfeval(get(my_handles{i}.plot_axes, 'Callback'), ...
%             my_handles{i}.plot_axes, []);
%     end
end
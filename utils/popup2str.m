function str = popup2str(popup_handle)
% str = popup2str(popup_handle)
%
% Returns the string currently selected in a popup box gui object.
%
c = cellstr(get(popup_handle, 'String'));
str = c{get(popup_handle, 'Value')};

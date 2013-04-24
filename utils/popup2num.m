function num = popup2num(popup_handle)
% num = popup2num(popup_handle)
%
% Returns the value of the integer selected in a popup box gui object.
%
c = cellstr(get(popup_handle, 'String'));
num = str2num(c{get(popup_handle, 'Value')});

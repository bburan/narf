function num = popup2num(handle)
% Returns the value of the integer selected in a popup box gui object.
    c = cellstr(get(handle, 'String'));
    num = str2num(c{get(handle, 'Value')});
end
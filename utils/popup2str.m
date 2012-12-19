function str = popup2str(handle)
% Returns the value of the integer selected in a popup box gui object.
    c = cellstr(get(handle, 'String'));
    str = c{get(handle, 'Value')};
end
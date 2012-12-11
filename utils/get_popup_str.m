function s = get_popup_str(module, guifield)
% Returns a string that is the value of the popup at that point.
    c = cellstr(get(module.plot_gui.(guifield), 'String'));
    s = c{get(module.plot_gui.(guifield), 'Value')};
end
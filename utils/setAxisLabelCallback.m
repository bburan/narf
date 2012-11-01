function [] = setAxisLabelCallback(gca, fn, axis)
% Sets up fn as the callback fn to label an axis dynamically
%
% INPUTS:
%   gca    Current graphics object handle
%   fn     Callback function which accepts a single argument only
%   axis   A string, either 'X' or 'Y' (case sensitive!)
%
% OUTPUTS:
%    none, since this is just used for its side effect
%
% 2012-09-21, Ivar Thorson

% Define an inner function because matlab is too retarded to make multiline anonymous ones
    function innerCallback(hProp,eventData)
        hAxes = eventData.AffectedObject;
        tickValues = get(hAxes, sprintf('%sTick', axis));
        newLabels = arrayfun(@(value)(fn(value)), tickValues, 'UniformOutput',false);
        set(hAxes, sprintf('%sTickLabel', axis), newLabels);
    end

hAxes = get(gcf,'CurrentAxes');
hhAxes = handle(hAxes);
hProp = findprop(hhAxes,  sprintf('%sTick', axis));
hListener = handle.listener(hhAxes, hProp, 'PropertyPostSet', @innerCallback);
setappdata(hAxes, sprintf('%sTickListener', axis)', hListener);
end

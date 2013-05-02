function [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack)
% [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack)
% 
% Returns the state of the plot controls from the first
% 'load_stim_resps_from_baphy' module. Used when determining what to plot.
% If a plot_gui cannot be found, it returns a safe default in case we want
% to use one of STACK's plot function from a script (ie, when there is
% no narf_modelpane gui connected so the model is headless).
%
% ARGUMENTS:
%    stack     A reference to the global STACK
%
% RETURNS:
%    sf        The selected stimfile to plot
%    stim_idx  The selected stimulus index number to plot
%    chan_idx  The selected channel index number to plot.
%

global XXX NARFGUI; 

[~, baphy_idx] = find_modules(stack, 'load_stim_resps_from_baphy', true);
if isempty(baphy_idx),
    [~,baphy_idx] = find_modules(stack, 'load_stim_resps_wehr',true);
end

if isfield(NARFGUI{baphy_idx}, 'plot_gui')
    sf = popup2str(NARFGUI{baphy_idx}.plot_gui.selected_stimfile_popup);
    stim_idx = popup2num(NARFGUI{baphy_idx}.plot_gui.selected_stim_idx_popup);
    chan_idx = popup2num(NARFGUI{baphy_idx}.plot_gui.selected_stim_chan_popup);
else
    sf = XXX{1}.training_set{1};
    stim_idx = 1;
    chan_idx = 1;
end


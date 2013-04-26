function delete_all_module_guis()
% delete_all_module_guis()
%
% Deletes any gui handles, widgets, etc from NARFGUI. This is important to do
% whenever the STACK is about to be destroyed, because otherwise you can
% lose access to the GUI handles and they will persist even when you want
% to create newer, replacement GUI widgets later.
%
% No arguments or return values. Use purely for side effects. 

global NARFGUI;

for ii = 1:length(NARFGUI)
    delete_module_gui(ii);
end

NARFGUI = {};
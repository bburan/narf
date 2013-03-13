function open_narf_gui()
 % Open a new NARF model pane figure to view the current model STACK 
 global XXX STACK NARF_MODULES_PATH MODULES;
 
 MODULES = scan_directory_for_modules(NARF_MODULES_PATH);

 pf = figure('Menubar','figure', 'Resize','off', 'MenuBar', 'none', ...
             'Units','pixels', 'Position', [20 50 1300 max(600, min(170*length(STACK)+40, 1100))]);
 narf_modelpane(pf); 
 
end
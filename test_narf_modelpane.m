% A function to test out narf_model.m 

w = 1300;  % Total pane width, including slider
h = 600;   % Total pane height

pf = figure('Menubar','figure', 'Resize','off', ...
    'Units','pixels', 'Position', [20 50 w h]);

narf_modelpane(pf, [],[]);

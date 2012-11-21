% A function to test out narf_model.m 

NARF_PATH       = '/home/ivar/matlab/narf/';
addpath([NARF_PATH filesep 'utils'], ...
        [NARF_PATH filesep 'modules']);


w = 1300;  % Total pane width, including slider
h = 600;   % Total pane height

pf = figure('Menubar','figure', 'Resize','off', ...
    'Units','pixels', 'Position', [20 50 w h]);

% Build a list of the files in the 'modules/' directory
mods = scan_directory_for_modules('~/matlab/narf/modules/');

% Pass modelpane a default starting STACK and X
narf_modelpane(pf, [], '', mods);
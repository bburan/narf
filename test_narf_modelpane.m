% A script to examine a single model by itself using 'narf_modelpane.m'
% If you aren't optimizing, this should be all you need.

NARF_PATH = '/home/ivar/matlab/narf/';
addpath([NARF_PATH filesep 'utils'], ...
        [NARF_PATH filesep 'modules']);

w = 1300;  % Total pane width, including slider
h = 600;   % Total pane height

pf = figure('Menubar','figure', 'Resize','off', ...
    'Units','pixels', 'Position', [20 50 w h]);

% Build a list of the files in the 'modules/' directory
mods = scan_directory_for_modules('~/matlab/narf/modules/');

stack = {};
stack{1} = [];
stack{1}.cellid = 'por012c-b1';
stack{1}.training_set = {'por012c02_p_TOR', 'por012c03_p_SPN'};
stack{1}.test_set = {};

% Pass modelpane the stack and X
narf_modelpane(pf, mods, stack, []); 

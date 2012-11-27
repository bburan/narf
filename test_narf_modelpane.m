% A script to examine a single model by itself using 'narf_modelpane.m'
% If you aren't optimizing, this should be all you need.

NARF_PATH = '/home/ivar/matlab/narf/';
addpath([NARF_PATH filesep 'utils'], ...
        [NARF_PATH filesep 'modules']);

pf = figure('Menubar','figure', 'Resize','off', ...
    'Units','pixels', 'Position', [20 50 1300 600]);

% Build a list of the files in the 'modules/' directory
mods = scan_directory_for_modules('~/matlab/narf/modules/');

global STACK XXX;

STACK = {};
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por012c-b1';
XXX{1}.training_set = {'por012c02_p_TOR', 'por012c03_p_SPN'};
XXX{1}.test_set = {};

% Pass modelpane the parent panel and all modules
narf_modelpane(pf, mods); 

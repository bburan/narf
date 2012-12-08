% function debugoff
%
% set global variable DEBUGON to 0... can be referenced by
% different .m files to switch from regular to debug mode
%
% SVD 11/8/01
function debugoff

global DEBUGON
DEBUGON=0;
disp('debugging is OFF.');

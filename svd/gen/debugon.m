% function debugon
%
% set global variable DEBUGON to 1... can be referenced by
% different .m files to switch from regular to debug mode
%
% SVD 11/8/01
function debugon

global DEBUGON
DEBUGON=1;
disp('debugging is ON.');

% function debugon=debugcheck()
%
% check the status of the global variable DEBUGON.  if not set,
% initialize to 0 (OFF)
%
% SVD 11/8/01
%
function debugon=debugcheck

global DEBUGON

if isempty(DEBUGON),
   DEBUGON=0;
end

debugon=DEBUGON;
return

if DEBUGON,
   disp('debugging is ON.');
else
   disp('debugging is OFF.');
end

debugon=DEBUGON;

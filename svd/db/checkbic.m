% function r=checkbic
%
% returns 1 if HOSTNAME contains ".bic.berkeley.edu"
% returns 2 if HOSTNAME contains ".Millenium.Berkeley.EDU"
% returns 0 otherwise
%
function r=checkbic

shost=getenv('HOSTNAME');
if findstr(shost,'.bic.berkeley.edu'),
   r=1;
elseif findstr(shost,'.Millennium.Berkeley.EDU'),
   r=2;
else
   r=0;
end


% function fullpage(orientation,fhandle)
%
function fullpage(orientation,fhandle)

if ~exist('orientation','var'),
   orientation='portrait';
end
if ~exist('fhandle','var'),
   fhandle=gcf;
end

if strcmp(lower(orientation),'portrait'),
   paperpos=[0.25 0.25 8 10.5];
elseif strcmp(lower(orientation),'landscape'),
   paperpos=[0.25 0.25 10.5 8];
else
   fprintf('valid orientation values are "portrait" or "landscape"');
end

set(fhandle,'PaperOrientation',orientation,'PaperPosition',paperpos);

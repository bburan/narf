% function y=hinge3(parms,x,oparms);
%
% parms=[x0 slope];
% x0 - threshold
% slope - global gain
%
function y=hinge3(parms,x,oparms);

if exist('oparms','var'),
   parms(find(~isnan(oparms)))=oparms(find(~isnan(oparms)));
end

phasecount=length(parms)-1;
spacecount=size(x,2);
x0=0;

if spacecount==1,
   
   slope=parms(2);
   y=(x-x0) .* slope;
elseif phasecount==1,
   
   slope=parms(2);
   y=(sum(x,2)-x0) .* slope;
elseif phasecount==spacecount,
   y=zeros(size(x,1),1);
   for phaseidx=1:phasecount,
      slope=parms(phaseidx+1);
      y=y+(x(:,phaseidx)-x0) .* slope;
   end
else
   y=zeros(size(x,1),1);
   stepsize=ceil(spacecount/phasecount);
   for phaseidx=1:phasecount,
      slope=parms(phaseidx+1);
      pidx=(1:stepsize)+stepsize*(phaseidx-1);
      pidx=pidx(find(pidx<=spacecount));
      tx=sum(x(:,pidx),2);
      y=y+(tx-x0) .* slope;
   end
end

% add dc after applying gain
x0=parms(1);
y=y+x0;

% apply threshold after summing over spatial channels
%y=y.*(y>0);

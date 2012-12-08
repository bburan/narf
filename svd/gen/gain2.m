% function y=gain2(parms,x,oparms);
%
% parms=[slope];
% slope - global gain
%
% ripped off hinge2 with x0 removed
%
function y=gain2(parms,x,oparms);

if exist('oparms','var'),
   parms(find(~isnan(oparms)))=oparms(find(~isnan(oparms)));
end

phasecount=length(parms);
spacecount=size(x,2);
if spacecount==1,
   slope=parms(1);
   
   y=x.*slope;
elseif phasecount==1,
   slope=parms(1);
   
   y=sum(x,2) .* slope;
elseif phasecount==spacecount,
   y=zeros(size(x,1),1);
   for phaseidx=1:phasecount,
      slope=parms(phaseidx);
      
      y=y+ x(:,phaseidx) .* slope;
   end
else
   y=zeros(size(x,1),1);
   stepsize=ceil(spacecount/phasecount);
   for phaseidx=1:phasecount,
      slope=parms(phaseidx);
      
      pidx=(1:stepsize)+stepsize*(phaseidx-1);
      pidx=pidx(find(pidx<=spacecount));
      tx=sum(x(:,pidx),2);
      y=y + tx .* slope;
   end
end


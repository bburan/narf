% function rt=thresh(parms,x);
%
% threshold a function.
%
% x     =  T x S vector (S=1 if pred has already been summed over
%          space, otherwise rectification happens in each spatial
%          channel and then gets summed)
% parms =  scalar threshold
%
% rt = sum( (x-parms) .* (x>parms) , 2 )
%
function rt=thresh(parms,x);

phasecount=length(parms);
spacecount=size(x,2);
if spacecount==1,
   rt=(x-parms(1)).*(x>=parms(1));
elseif phasecount==1,
   rt=sum((x-parms).*(x>parms),2);
else
   rt=zeros(size(x,1),1);
   stepsize=ceil(spacecount/phasecount);
   for phaseidx=1:phasecount,
      pidx=(1:stepsize)+stepsize*(phaseidx-1);
      pidx=pidx(find(pidx<=spacecount));
      tx=sum(x(:,pidx),2);
      rt=rt+(tx-parms(phaseidx)).*(tx>=parms(phaseidx));
   end
end
